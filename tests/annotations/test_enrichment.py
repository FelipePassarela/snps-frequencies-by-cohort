from pathlib import Path

import pandas as pd

from usp_gene_correlation.annotations import CsvAnnotationSource, annotate_frequencies
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
