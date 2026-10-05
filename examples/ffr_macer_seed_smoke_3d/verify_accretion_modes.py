#!/usr/bin/env python3
import argparse
import math
import re
from pathlib import Path

HERE = Path(__file__).resolve().parent
CASES = {
    "tng_bondi": (0, "tng-bondi"),
    "boosted_bondi": (1, "boosted-bondi"),
    "am_bondi": (2, "am-bondi"),
    "ffr": (3, "ffr"),
}
KV = re.compile(r"([A-Za-z_][A-Za-z0-9_]*)=([^\s]+)")

def fields(line):
    return dict(KV.findall(line))

def fval(d, key):
    try:
        x = float(d[key])
    except Exception as exc:
        raise AssertionError(f"missing/invalid {key} in {d}") from exc
    if not math.isfinite(x):
        raise AssertionError(f"non-finite {key}={x}")
    return x

def parse_param(path):
    out = {}
    for raw in path.read_text().splitlines():
        s = raw.strip()
        if not s or s.startswith("%"):
            continue
        parts=s.split()
        if len(parts) >= 2:
            out[parts[0]]=parts[1]
    return out

def relerr(a,b):
    return abs(a-b)/max(abs(a),abs(b),1e-300)

def check_case(name, model, model_name, args):
    log = HERE / f"run_acc_{name}.log"
    param = HERE / f"param_acc_{name}.txt"
    if not log.is_file() or not param.is_file():
        raise AssertionError(f"{name}: missing log or parameter file")

    p = parse_param(param)
    assert int(p["BHBenchmarkAccretionModel"]) == model, (name, p["BHBenchmarkAccretionModel"])
    assert int(p["BHBenchmarkAccretionTarget"]) == 1, (name, p["BHBenchmarkAccretionTarget"])
    assert int(p["BHBenchmarkFeedbackModel"]) == 0, (name, p["BHBenchmarkFeedbackModel"])
    assert math.isclose(float(p["BHBenchmarkEddingtonFactor"]), args.expected_edd_factor,
                        rel_tol=1e-12, abs_tol=0.0), (name, p["BHBenchmarkEddingtonFactor"])

    text = log.read_text(errors="replace")
    lowered = text.lower()
    for bad in ("terminate(", "termination condition", "fatal error", "ledger failed"):
        if bad in lowered:
            raise AssertionError(f"{name}: failure marker in log: {bad}")
    if "BH_TNG: thermal event" in text or "BH_TNG: kinetic event" in text or "BH_TNG: source" in text:
        raise AssertionError(f"{name}: TNG feedback activity present with feedback=NONE")
    if "BH_FFR: wind packet fired" in text or "BH_FFR: jet packet fired" in text:
        raise AssertionError(f"{name}: MACER packet fired in direct+NONE run")

    acc = [fields(x) for x in text.splitlines() if "BH_BENCHMARK: accretion" in x]
    if not (args.min_transactions <= len(acc) <= args.max_transactions):
        raise AssertionError(f"{name}: transaction count {len(acc)} outside "
                             f"[{args.min_transactions},{args.max_transactions}]")

    raws=[]
    for d in acc:
        if d.get("model") != model_name:
            raise AssertionError(f"{name}: wrong model in log: {d.get('model')}")
        if d.get("target") != "direct":
            raise AssertionError(f"{name}: wrong target in log: {d.get('target')}")
        if int(d.get("feedback","-1")) != 0:
            raise AssertionError(f"{name}: wrong feedback selector in log: {d.get('feedback')}")
        raw=fval(d,"raw"); op=fval(d,"operational")
        if raw < 0 or op < 0:
            raise AssertionError(f"{name}: negative raw/operational rate")
        if raw > 0 and relerr(raw,op) > args.raw_operational_tol:
            raise AssertionError(f"{name}: estimator-only cap is active: raw={raw} op={op}")
        raws.append(raw)

    positive=[x for x in raws if x > 0]
    if not positive:
        raise AssertionError(f"{name}: no positive raw estimator samples")
    med=sorted(positive)[len(positive)//2]
    drift=(max(positive)-min(positive))/max(abs(med),1e-300)
    if drift > args.raw_stability_tol:
        raise AssertionError(f"{name}: raw-rate short-run drift {drift:.6g} > {args.raw_stability_tol}")

    sink=[fields(x) for x in text.splitlines() if "BH_BENCHMARK: sink" in x]
    physical=[]
    full_active=[]
    hetero=0
    partial=0
    for d in sink:
        r=fval(d,"R_sink")
        active_frac=fval(d,"activeGasFrac")
        cells=int(d.get("sinkCells","0"))
        if int(d.get("heterogeneous","0")):
            hetero += 1
        if cells > 0:
            physical.append(r)
            if active_frac >= args.full_active_fraction:
                full_active.append(r)
            else:
                partial += 1
    if not physical:
        raise AssertionError(f"{name}: no physical sink transaction with sinkCells>0")
    if not full_active:
        raise AssertionError(f"{name}: no transaction had a fully active accretion aperture; "
                             "cannot apply the per-transaction R_sink unity check")
    worst=max(abs(1.0-r) for r in full_active)
    if worst > args.rsink_tol:
        raise AssertionError(f"{name}: fully-active worst |1-R_sink|={worst:.6g} > {args.rsink_tol}")

    ledgers=[fields(x) for x in text.splitlines() if "BH_BENCHMARK: mass-ledger" in x]
    if not ledgers:
        raise AssertionError(f"{name}: no cumulative mass ledger")
    last=ledgers[-1]
    cop=fval(last,"cumOperational")
    creal=fval(last,"cumRealized")
    cbh=fval(last,"cumBHGrowth")
    if relerr(creal,cbh) > args.mass_closure_tol:
        raise AssertionError(f"{name}: realized/BH-growth closure failed: {creal} vs {cbh}")
    cumulative_mismatch = relerr(cop,creal) if cop > 0 or creal > 0 else 0.0

    print(f"{name:16s} N={len(acc):2d} raw_drift={drift:.3e} "
          f"full-active max|1-R_sink|={worst:.3e} cumulative_op_vs_real={cumulative_mismatch:.3e} "
          f"heterogeneous_sink_txn={hetero} partial-active_txn={partial}")
    return {
        "transactions": len(acc),
        "raw_drift": drift,
        "worst_rsink": worst,
        "cumulative_mismatch": cumulative_mismatch,
        "heterogeneous": hetero,
        "partial_active": partial,
    }

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--rsink-tol", type=float, default=1e-3)
    ap.add_argument("--raw-stability-tol", type=float, default=5e-2)
    ap.add_argument("--raw-operational-tol", type=float, default=1e-10)
    ap.add_argument("--mass-closure-tol", type=float, default=2e-10)
    ap.add_argument("--expected-edd-factor", type=float, default=1e30)
    ap.add_argument("--full-active-fraction", type=float, default=1.0-1e-8)
    ap.add_argument("--min-transactions", type=int, default=2)
    ap.add_argument("--max-transactions", type=int, default=12)
    args=ap.parse_args()

    for name,(model,model_name) in CASES.items():
        check_case(name,model,model_name,args)
    print("Accretion-only benchmark verification: PASS")

if __name__ == "__main__":
    main()