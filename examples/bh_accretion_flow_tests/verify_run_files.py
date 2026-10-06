#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

HERE = Path(__file__).resolve().parent

REQUIRED_CONFIG = {
    "DOUBLEPRECISION=1",
    "DOUBLEPRECISION_FFTW",
    "INPUT_IN_DOUBLEPRECISION",
    "OUTPUT_IN_DOUBLEPRECISION",
    "OUTPUT_COORDINATES_IN_DOUBLEPRECISION",
    "NGB_TREE_DOUBLEPRECISION",
    "HAVE_HDF5",
    "BLACKHOLE_FFR",
    "FOF",
}

REQUIRED_PARAMS = {
    "ICFormat": "3",
    "SnapFormat": "3",
    "ComovingIntegrationOn": "0",
    "PeriodicBoundariesOn": "1",
    "BoxSize": "0.1",
    "OutputListOn": "1",
    "BHAccretionRadius": "0.0125",
    "BHFeedbackRadius": "0.025",
    "BHFreeFallA": "1.0e-3",
    "BHBenchmarkAccretionModel": "4",
    "BHBenchmarkAccretionTarget": "1",
    "BHBenchmarkFeedbackModel": "0",
    "BHBenchmarkEddingtonFactor": "1.0e-8",
}


def parse_param(path: Path) -> dict[str, str]:
    out = {}
    for raw in path.read_text().splitlines():
        s = raw.strip()
        if not s or s.startswith("%"):
            continue
        parts = s.split()
        if len(parts) >= 2:
            out[parts[0]] = parts[1]
    return out


def main() -> None:
    cfg = HERE / "Config_flow_tests.sh"
    config_tokens = {line.strip() for line in cfg.read_text().splitlines()
                     if line.strip() and not line.lstrip().startswith("#")}
    missing = REQUIRED_CONFIG - config_tokens
    if missing:
        raise SystemExit("Config missing: " + ", ".join(sorted(missing)))

    manifest_path = HERE / "run_manifest.json"
    if not manifest_path.exists():
        raise SystemExit("Missing run_manifest.json; run python3 generate_run_files.py")
    manifest = json.loads(manifest_path.read_text())

    if manifest["diagnostic_count"] != 42:
        raise SystemExit(f"Expected 42 diagnostic runs, got {manifest['diagnostic_count']}")
    if manifest["evolved_count"] != 18:
        raise SystemExit(f"Expected 18 evolved runs, got {manifest['evolved_count']}")
    if manifest["total_count"] != 60:
        raise SystemExit(f"Expected 60 total planned runs, got {manifest['total_count']}")

    outputlists = set()
    for run in manifest["runs"]:
        p = HERE / run["parameter_file"]
        if not p.exists():
            raise SystemExit(f"Missing parameter file {p}")
        d = parse_param(p)
        for key, expected in REQUIRED_PARAMS.items():
            if d.get(key) != expected:
                raise SystemExit(f"{p}: {key}={d.get(key)!r}, expected {expected!r}")

        ic = HERE / (d["InitCondFile"] + ".hdf5")
        if not ic.exists():
            raise SystemExit(f"{p}: missing IC {ic}")

        ol = HERE / d["OutputListFilename"]
        if not ol.exists():
            raise SystemExit(f"{p}: missing output list {ol}")
        outputlists.add(ol)

        times = []
        for line in ol.read_text().splitlines():
            if not line.strip():
                continue
            fields = line.split()
            if len(fields) != 2 or fields[1] != "2":
                raise SystemExit(f"{ol}: every output must use DumpFlag=2, got {line!r}")
            times.append(float(fields[0]))
        if not times:
            raise SystemExit(f"{ol}: empty output list")
        tmax = float(d["TimeMax"])
        if abs(times[-1] - tmax) > 5e-13 * max(abs(tmax), 1e-30):
            raise SystemExit(f"{p}: final output {times[-1]} != TimeMax {tmax}")

    print("Run-file verification: PASS")
    print("  diagnostics: 42")
    print("  representative evolved: 18")
    print("  total parameter files: 60")
    print(f"  shared output lists: {len(outputlists)}")
    print("  all outputs use DumpFlag=2 (full snapshot, no FoF catalogue)")


if __name__ == "__main__":
    main()
