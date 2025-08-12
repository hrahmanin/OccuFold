import argparse, os, re, sys
import pandas as pd

# Add the repo root (parent of this scripts/ folder) to sys.path
REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir))
if REPO_ROOT not in sys.path:
    sys.path.insert(0, REPO_ROOT)

# Try your local package first, then a flat utils/ as fallback
def import_modules():
    try:
        import OccuFold.utils.snippet as snp
        import OccuFold.utils.experimental_path as exp
        return snp, exp, "OccuFold.utils"
    except ModuleNotFoundError:
        try:
            import utils.snippet as snp
            import utils.experimental_path as exp
            return snp, exp, "utils"
        except ModuleNotFoundError as e:
            print("Could not import your utils module. sys.path:", file=sys.stderr)
            for p in sys.path:
                print(" -", p, file=sys.stderr)
            raise e

def sanitize(name: str) -> str:
    return re.sub(r'[^A-Za-z0-9._-]', '_', name)

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--region", required=True)
    p.add_argument("--outdir", required=True)
    p.add_argument("--peaks", default=None)
    p.add_argument("--motifs", default=None)
    a = p.parse_args()

    os.makedirs(a.outdir, exist_ok=True)
    os.makedirs("processing", exist_ok=True)

    snp, exp, used = import_modules()
    print(f"[step1] Using imports from: {used}")

    peaks  = a.peaks  if a.peaks  else exp.ctcf_peaks
    motifs = a.motifs if a.motifs else exp.ctcf_motifs

    df = snp.get_region_snippet(peaks, motifs, a.region)
    df = pd.DataFrame(df)

    base = sanitize(a.region)
    p1 = os.path.join("processing", f"{base}.csv")
    p2 = os.path.join(a.outdir,     f"{base}.csv")
    df.to_csv(p1, index=False)
    df.to_csv(p2, index=False)
    print(f"Wrote {p1} and {p2}")

if __name__ == "__main__":
    main()
