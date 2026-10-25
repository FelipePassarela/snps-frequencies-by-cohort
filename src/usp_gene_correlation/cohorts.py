"""Read cohort sample identifiers."""
import re
from pathlib import Path
import pandas as pd

def get_cohort_ids(cohorts_path: list[Path]) -> pd.DataFrame:
    """Get cohort IDs from CSV files.

    Args:
        cohorts_path (list[Path]): List of paths to cohort CSV files.

    Returns:
        pd.DataFrame: DataFrame with cohort IDs as columns.

    Raises:
        ValueError: If no matching columns are found in a cohort.
    """

    if not cohorts_path:
        raise ValueError("At least one cohort is required")
    if len({path.stem for path in cohorts_path}) != len(cohorts_path):
        raise ValueError("Cohort filenames must have unique stems")
    cohorts: list[pd.DataFrame] = [
        pd.read_csv(path, sep=None, engine="python", header=None)
        for path in cohorts_path
    ]
    max_len = max(len(cohort) for cohort in cohorts)

    for i, cohort in enumerate(cohorts):
        # The files do not have a consistent header, so we search for the ID column
        id_pattern = re.compile(r"^C\d+-ExC\d+-xgenV\d+$")
        col_matches = [
            col
            for col in cohort.columns
            if cohort[col].astype(str).str.match(id_pattern).any()
        ]
        id_col = next(iter(col_matches), None)
        if id_col is None:
            raise ValueError(
                f"No matching {id_pattern.pattern} columns found in cohort."
            )

        # Only ID column is required
        cohort = cohort.loc[cohort[id_col].astype(str).str.match(id_pattern), [id_col]].drop_duplicates().reset_index(drop=True)
        cohorts[i] = cohort

        # Some cohorts have fewer IDs; pad with NaNs
        if len(cohort) < max_len:
            diff = max_len - len(cohort)
            padding = pd.DataFrame({id_col: [pd.NA] * diff})
            cohorts[i] = pd.concat([cohort, padding], ignore_index=True)

    cohort_ids = pd.concat(cohorts, axis=1)
    cohort_ids.columns = [path.stem for path in cohorts_path]
    return cohort_ids


