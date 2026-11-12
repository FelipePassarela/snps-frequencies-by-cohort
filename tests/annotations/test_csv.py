from pathlib import Path

import pytest

from usp_gene_correlation.annotations import CsvAnnotationSource, SnpAnnotation


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
