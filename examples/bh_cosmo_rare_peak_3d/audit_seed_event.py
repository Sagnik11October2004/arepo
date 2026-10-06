#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

SEED_RE = re.compile(
    r"BH_FFR: seeded FoF group (?P<group>\d+) at z=(?P<z>[-+0-9.eE]+) "
    r"with Mhalo=(?P<mhalo>[-+0-9.eE]+) Msun, MBH=(?P<mbh>[-+0-9.eE]+) Msun, "
    r"central gas ID=(?P<bhid>\d+)\."
)
SELECT_RE = re.compile(
    r"BH_FFR: single-seed mode selected most massive eligible FoF group "
    r"(?P<group>\d+) with Mhalo=(?P<mhalo>[-+0-9.eE]+) Msun\."
)
COUNT_RE = re.compile(
    r"BH_FFR: FoF full-step seeding created (?P<created>\d+) seed\(s\); "
    r"total Type-5 BH count is now (?P<total>\d+)\."
)


def main() -> None:
    ap = argparse.ArgumentParser(
        description="Audit the controlled rare-peak FoF seed transaction."
    )
    ap.add_argument("log", type=Path)
    ap.add_argument("output", type=Path)
    ap.add_argument("--z-min", type=float, default=20.0)
    ap.add_argument("--z-max", type=float, default=22.0)
    ap.add_argument("--expected-seed-mass", type=float, default=1.0e5)
    args = ap.parse_args()

    text = args.log.read_text(errors="replace")
    seeds = [m.groupdict() for m in SEED_RE.finditer(text)]
    if len(seeds) != 1:
        raise SystemExit(
            f"Expected exactly one FoF seed event in {args.log}; found {len(seeds)}"
        )

    seed = seeds[0]
    group = int(seed["group"])
    z = float(seed["z"])
    mhalo = float(seed["mhalo"])
    mbh = float(seed["mbh"])
    bhid = int(seed["bhid"])

    if not (args.z_min < z < args.z_max):
        raise SystemExit(
            f"Seed redshift is outside the strict requested window: "
            f"z={z:.9g}, require {args.z_min}<z<{args.z_max}"
        )
    if abs(mbh / args.expected_seed_mass - 1.0) > 2.0e-6:
        raise SystemExit(
            f"Unexpected seed mass in log: {mbh:.9g} Msun; "
            f"expected {args.expected_seed_mass:.9g} Msun"
        )

    selected = [m.groupdict() for m in SELECT_RE.finditer(text)]
    if selected:
        matching = [
            s for s in selected
            if int(s["group"]) == group
            and abs(float(s["mhalo"]) / mhalo - 1.0) < 2.0e-6
        ]
        if not matching:
            raise SystemExit(
                "Seed event does not match any logged single-most-massive selection."
            )

    counts = [m.groupdict() for m in COUNT_RE.finditer(text)]
    one_seed_counts = [
        c for c in counts if int(c["created"]) == 1 and int(c["total"]) == 1
    ]
    if not one_seed_counts:
        raise SystemExit(
            "Did not find a FoF completion record confirming exactly one live Type-5 BH."
        )

    record = {
        "bh_id": bhid,
        "fof_group": group,
        "host_fof_mass_msun": mhalo,
        "seed_mass_msun": mbh,
        "seed_redshift": z,
        "seed_window": {
            "strict_min_redshift": args.z_min,
            "strict_max_redshift": args.z_max,
        },
        "source_log": str(args.log),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(record, indent=2, sort_keys=True) + "\n")

    print(
        "Validated seed event: "
        f"z={z:.6f}, FoF={group}, Mhalo={mhalo:.6e} Msun, "
        f"MBH={mbh:.6e} Msun, BH_ID={bhid}"
    )
    print(f"Wrote seed metadata: {args.output}")


if __name__ == "__main__":
    main()
