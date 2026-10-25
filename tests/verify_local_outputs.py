"""Manual regression check against local (untracked) study data."""
from pathlib import Path
from tempfile import TemporaryDirectory

import pandas as pd
from usp_gene_correlation.cli import main

def verify_and_update() -> None:
    cohorts = sorted(Path("data/cohorts").glob("*.csv"))
    verified: list[tuple[Path, bytes]] = []
    with TemporaryDirectory() as temporary:
        for group in ("longa", "grave"):
            out = Path(temporary) / group
            main([f"data/{group}.vcf", *(str(p) for p in cohorts), "-o", str(out)])
            for generated in sorted(out.glob("*.csv")):
                target = Path("output") / group / generated.name
                before = pd.read_csv(target, index_col=0).drop(columns=["name", "rsid"], errors="ignore")
                after = pd.read_csv(generated, index_col=0)
                pd.testing.assert_frame_equal(before, after.drop(columns=["name", "rsid"]))
                assert after[["name", "rsid"]].notna().all().all()
                verified.append((target, generated.read_bytes()))
        for target, content in verified:
            target.write_bytes(content)
    print(f"Updated {len(verified)} CSVs; all original frequencies preserved.")

if __name__ == "__main__":
    verify_and_update()
