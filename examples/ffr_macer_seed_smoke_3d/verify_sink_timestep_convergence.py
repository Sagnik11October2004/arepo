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

def metrics(log):
    text=log.read_text(errors="replace")
    sink=[fields(x) for x in text.splitlines() if "BH_BENCHMARK: sink" in x]
    physical=[]
    hetero=0
    for d in sink:
        cells=int(d.get("sinkCells","0"))
        if cells > 0:
            r=float(d["R_sink"])
            if not math.isfinite(r):
                raise AssertionError(f"{log}: non-finite R_sink")
            physical.append(r)
            hetero += int(d.get("heterogeneous","0") != "0")
    if not physical:
        raise AssertionError(f"{log}: no physical sink samples")

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
        "maxerr": max(abs(1-r) for r in physical),
        "meanerr": sum(abs(1-r) for r in physical)/len(physical),
        "cumerr": relerr(cop,creal) if cop > 0 or creal > 0 else 0.0,
        "hetero": hetero,
        "n": len(physical),
    }

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--root", required=True)
    ap.add_argument("--step-fracs", nargs="+", type=float, required=True)
    ap.add_argument("--finest-rsink-tol", type=float, default=1e-3)
    ap.add_argument("--trend-slack", type=float, default=0.10)
    args=ap.parse_args()

    root=Path(args.root)
    fracs=sorted(args.step_fracs, reverse=True)
    allm={case:[] for case in CASES}

    print("model            step_frac    N   max|1-R|    mean|1-R|   cumulative_mass_err  hetero")
    for frac in fracs:
        d=root/f"step_{tag(frac)}"
        for case in CASES:
            m=metrics(d/f"run_acc_{case}.log")
            allm[case].append((frac,m))
            print(f"{case:16s} {frac:10.6g} {m['n']:3d} {m['maxerr']:11.3e} "
                  f"{m['meanerr']:11.3e} {m['cumerr']:19.3e} {m['hetero']:7d}")

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
