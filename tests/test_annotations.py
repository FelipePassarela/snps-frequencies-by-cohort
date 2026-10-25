from pathlib import Path

import pandas as pd
import pytest

from usp_gene_correlation.annotations import (
    CsvAnnotationSource,
    SnpAnnotation,
    annotate_frequencies,
)
from usp_gene_correlation.genotypes import GenotypeFrequencies


def test_annotations_preserve_values_order_and_unknowns(tmp_path: Path) -> None:
    original = pd.DataFrame(
        {
            "0/0": [0.25, 1.0, 0.5],
            "0/1": [0.75, 0, 0.5],
            "1/1": [0, 0, 0],
            "./.": [0, 0, 0],
        },
        index=["chr2_2", "1_1", "chr2_2"],
    )
    frequencies = GenotypeFrequencies(original)
    catalog = tmp_path / "annotations.csv"
    catalog.write_text("chr_pos,name,rsid\n2_2,GENE,rs1\n", encoding="utf-8")
    source = CsvAnnotationSource(catalog)
    actual = annotate_frequencies(frequencies, source)
    pd.testing.assert_frame_equal(
        actual.drop(columns=["name", "rsid"]), frequencies.to_dataframe()
    )
    assert actual.name.tolist() == ["GENE", "", "GENE"]
    assert actual.rsid.tolist() == ["rs1", "", "rs1"]
    assert list(original.columns) == ["0/0", "0/1", "1/1", "./."]


@pytest.mark.parametrize(
    "contents",
    [
        "position,name,rsid\nchr1_1,A,rs1\n",
        "chr_pos,name,rsid\nchr1_1,A,rs1\n1_1,B,rs2\n",
        "chr_pos,name,rsid\nchr1_0,A,rs1\n",
        "chr_pos,name,rsid\nchr1_1,A,bad\n",
        "chr_pos,name,rsid\nchr1_1,A\n",
    ],
)
def test_invalid_catalog(tmp_path: Path, contents: str) -> None:
    path = tmp_path / "snps_annotations.csv"
    path.write_text(contents)
    with pytest.raises(ValueError):
        CsvAnnotationSource(path)


def test_bundled_catalog(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.chdir(tmp_path)
    source = CsvAnnotationSource.from_default_catalog()
    assert source.lookup("9_62897604") == SnpAnnotation(
        "chr9_62897604", "LERFS", "rs200443822"
    )
    assert source.lookup("chr2_46818170").name == "LINC01118"
    assert len(source._records) == 44
