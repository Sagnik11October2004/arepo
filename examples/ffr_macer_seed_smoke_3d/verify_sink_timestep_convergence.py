#!/usr/bin/env python3
import argparse
import math
import re
from pathlib import Path

CASES=("tng_bondi","boosted_bondi","am_bondi","ffr")
KV=re.compile(r"([A-Za-z_][A-Za-z0-9_]*)=([^\s]+)")

def fields(line):
    return dict(KV.findall(line))

def relerr(a,b):
    return abs(a-b)/max(abs(a),abs(b),1e-300)

def tag(frac):
    return ("%g" % frac).replace(".","p")

def metrics(log, full_active_fraction):
    text=log.read_text(errors="replace")
    sink=[fields(x) for x in text.splitlines() if "BH_BENCHMARK: sink" in x]
    physical=[]
    full_active=[]
    hetero=0
    partial=0
    for d in sink:
        cells=int(d.get("sinkCells","0"))
        if cells > 0:
            r=float(d["R_sink"])
            af=float(d["activeGasFrac"])
            if not math.isfinite(r) or not math.isfinite(af):
                raise AssertionError(f"{log}: non-finite sink diagnostic")
            physical.append(r)
            if af >= full_active_fraction:
                full_active.append(r)
            else:
                partial += 1
            hetero += int(d.get("heterogeneous","0") != "0")
    if not physical:
        raise AssertionError(f"{log}: no physical sink samples")
    if not full_active:
        raise AssertionError(f"{log}: no fully active aperture transaction for instantaneous R_sink test")

    led=[fields(x) for x in text.splitlines() if "BH_BENCHMARK: mass-ledger" in x]
    if not led:
        raise AssertionError(f"{log}: no mass ledger")
    last=led[-1]
    cop=float(last["cumOperational"])
    creal=float(last["cumRealized"])
    cbh=float(last["cumBHGrowth"])
    if relerr(creal,cbh) > 2e-10:
        raise AssertionError(f"{log}: direct BH-growth ledger does not close")

    return {
        "maxerr": max(abs(1-r) for r in full_active),
        "meanerr": sum(abs(1-r) for r in full_active)/len(full_active),
        "cumerr": relerr(cop,creal) if cop > 0 or creal > 0 else 0.0,
        "hetero": hetero,
        "partial": partial,
        "n": len(full_active),
    }

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--root", required=True)
    ap.add_argument("--step-fracs", nargs="+", type=float, required=True)
    ap.add_argument("--finest-rsink-tol", type=float, default=1e-3)
    ap.add_argument("--trend-slack", type=float, default=0.10)
    ap.add_argument("--full-active-fraction", type=float, default=1.0-1e-8)
    args=ap.parse_args()

    root=Path(args.root)
    fracs=sorted(args.step_fracs, reverse=True)
    allm={case:[] for case in CASES}

    print("model            step_frac  Nfull max|1-R|    mean|1-R|   cumulative_mass_err  hetero partial")
    for frac in fracs:
        d=root/f"step_{tag(frac)}"
        for case in CASES:
            m=metrics(d/f"run_acc_{case}.log", args.full_active_fraction)
            allm[case].append((frac,m))
            print(f"{case:16s} {frac:10.6g} {m['n']:3d} {m['maxerr']:11.3e} "
                  f"{m['meanerr']:11.3e} {m['cumerr']:19.3e} {m['hetero']:7d} {m['partial']:7d}")

    for case, rows in allm.items():
        coarse=rows[0][1]
        fine=rows[-1][1]
        if fine["maxerr"] > args.finest_rsink_tol:
            raise AssertionError(f"{case}: finest max |1-R_sink|={fine['maxerr']} > {args.finest_rsink_tol}")
        # The exact exponential sink should approach the instantaneous hazard
        # as timesteps shrink. Permit small state/roundoff noise rather than
        # demanding strict monotonicity at every intermediate level.
        if fine["maxerr"] > coarse["maxerr"]*(1+args.trend_slack) + 1e-12:
            raise AssertionError(f"{case}: R_sink did not improve from coarse to fine")
        if fine["cumerr"] > coarse["cumerr"]*(1+args.trend_slack) + 1e-12:
            raise AssertionError(f"{case}: cumulative operational/realized mismatch did not improve")

    hetero_total=sum(m["hetero"] for rows in allm.values() for _,m in rows)
    if hetero_total:
        print(f"Heterogeneous hydro-timebin sink transactions observed: {hetero_total}")
    else:
        print("No heterogeneous hydro-timebin sink transaction occurred in this source; "
              "the dedicated L0/L1/L2 ICs should exercise this separately.")
    print("Sink timestep convergence verification: PASS")

if __name__ == "__main__":
    main()