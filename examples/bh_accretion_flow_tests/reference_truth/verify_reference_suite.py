#!/usr/bin/env python3
"""Static verification for generated high-resolution reference ICs and params."""
from __future__ import annotations

import json
from pathlib import Path

import h5py
import numpy as np


HERE = Path(__file__).resolve().parent
FLOW = HERE.parent


def require(cond: bool, message: str) -> None:
    if not cond:
        raise AssertionError(message)


def parse_params(path: Path) -> dict[str, str]:
    out = {}
    for raw in path.read_text().splitlines():
        line = raw.strip()
        if not line or line.startswith("%"):
            continue
        parts = line.split()
        if len(parts) >= 2:
            out[parts[0]] = parts[1]
    return out


def main() -> None:
    icm = json.loads((HERE / "reference_ic_manifest.json").read_text())["ics"]
    runm = json.loads((HERE / "reference_run_manifest.json").read_text())["runs"]

    require(len(icm) == 42, f"expected 42 reference ICs, got {len(icm)}")
    require(len(runm) == 42, f"expected 42 reference runs, got {len(runm)}")

    expected_per_level = {"HR4": 14, "HR5": 14, "HR6": 14}
    for level, expected in expected_per_level.items():
        n = sum(x["level"] == level for x in icm)
        require(n == expected, f"{level}: expected {expected} ICs, got {n}")

    for x in icm:
        p = FLOW / x["ic"]
        require(p.exists(), f"missing IC {p}")
        require(abs(x["recommended_sink_radius_pc"] - 4.0 * x["finest_dx_pc"]) < 1e-12,
                f"{p}: sink radius is not 4 finest cells")

        with h5py.File(p, "r") as f:
            n0 = int(f["Header"].attrs["NumPart_Total"][0])
            require(n0 == x["ngas"], f"{p}: gas count mismatch")
            for name in ("Coordinates", "Velocities", "Masses", "InternalEnergy"):
                require(f["PartType0"][name].dtype == np.dtype("float64"),
                        f"{p}: {name} is not float64")
            require(f["PartType5"]["Coordinates"].dtype == np.dtype("float64"),
                    f"{p}: BH coordinates are not float64")

    for r in runm:
        p = FLOW / r["parameter_file"]
        require(p.exists(), f"missing params {p}")
        prm = parse_params(p)

        require(prm["BHBenchmarkAccretionModel"] == "3", f"{p}: absorber must use volume FFR backend")
        require(prm["BHBenchmarkAccretionTarget"] == "1", f"{p}: absorber must be direct")
        require(prm["BHBenchmarkFeedbackModel"] == "0", f"{p}: feedback must be off")
        require(float(prm["BHBenchmarkEddingtonFactor"]) >= 1e20, f"{p}: Eddington cap is not passive")
        require(abs(float(prm["BHFreeFallA"]) - float(r["absorber_A"])) < 1e-12,
                f"{p}: absorber A mismatch")
        require(abs(float(prm["BHMaxSinkFraction"]) - float(r["absorber_max_fraction"])) < 1e-12,
                f"{p}: sink fraction mismatch")
        require(abs(1000.0 * float(prm["BHAccretionRadius"]) - float(r["sink_radius_pc"])) < 1e-10,
                f"{p}: sink radius mismatch")
        require(float(prm["SofteningComovingType1"]) <= 2.0e-5,
                f"{p}: BH softening too large for reference sink")

        output_list = FLOW / r["output_list"]
        require(output_list.exists(), f"missing output list {output_list}")
        lines = [q.split() for q in output_list.read_text().splitlines() if q.strip()]
        require(len(lines) == 20, f"{output_list}: expected 20 output times")
        require(all(len(q) >= 2 and q[1] == "2" for q in lines),
                f"{output_list}: every output must use DumpFlag=2")

    print("Reference-suite verification: PASS")
    print("  ICs: 42")
    print("  runs: 42")
    print("  levels: HR4/HR5/HR6 = 14 each")
    print("  sink radius: 4 finest cells")
    print("  absorber: direct, no feedback, passive Eddington cap")
    print("  outputs: 20 snapshots/run, DumpFlag=2")


if __name__ == "__main__":
    main()
