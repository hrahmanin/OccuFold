# scripts/step3_barriers.py
#!/usr/bin/env python3
import argparse, os, json
import pandas as pd
import numpy as np

from OccuFold.utils import convert

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--input", required=True, help="CSV from step2 (*.occupancy.csv)")
    p.add_argument("--region", required=True)
    p.add_argument("--outdir", required=True)
    p.add_argument("--lattice-site", type=int, default=250)
    p.add_argument("--paramdict", default=None, help="Optional JSON with simulation params")
    args = p.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    df = pd.read_csv(args.input)
    paramdict = {}
    if args.paramdict and os.path.exists(args.paramdict):
        with open(args.paramdict) as fh:
            paramdict = json.load(fh)
    # Compute and store lattice-derived settings
    lat_size = int(convert.get_lattice_size(args.region, lattice_site=args.lattice_site))
    paramdict["lattice_site"] = args.lattice_site            # provenance
    paramdict["region"] = args.region                        # provenance
    paramdict["monomers_per_replica"] = lat_size // 10       # your rule

    refined_occupancy = convert.get_refined_occupancy(df, args.region)

    CTCF_left_positions, CTCF_right_positions, ctcf_loc_list, ctcf_lifetime_list, ctcf_offtime_list = convert.get_ctcf_list(
        refined_occupancy, paramdict
    )

    paramdict['monomers_per_replica'] = int(
        convert.get_lattice_size(args.region, lattice_site=args.lattice_site) // 10
    )

    stem = os.path.splitext(os.path.basename(args.input))[0]

    ro_df = refined_occupancy if isinstance(refined_occupancy, pd.DataFrame) else pd.DataFrame(refined_occupancy)
    ro_df.to_csv(os.path.join(args.outdir, f"{stem}.refined_occupancy.csv"), index=False)

    barriers = pd.DataFrame({
        "side": (["L"] * len(CTCF_left_positions)) + (["R"] * len(CTCF_right_positions)),
        "position": list(map(int, CTCF_left_positions)) + list(map(int, CTCF_right_positions))
    })
    barriers.to_csv(os.path.join(args.outdir, f"{stem}.barriers.csv"), index=False)

    lists = pd.DataFrame({
        "loc":       list(map(int, ctcf_loc_list)),
        "lifetime":  list(map(float, ctcf_lifetime_list)),
        "offtime":   list(map(float, ctcf_offtime_list)),
    })
    lists.to_csv(os.path.join(args.outdir, f"{stem}.ctcf_lists.csv"), index=False)

    with open(os.path.join(args.outdir, f"{stem}.paramdict.json"), "w") as fh:
        json.dump(paramdict, fh, indent=2)

if __name__ == "__main__":
    main()
