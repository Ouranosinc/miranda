from pathlib import Path

import pandas as pd
import pytest  # noqa

from miranda.structure import build_path_from_schema


def get_testing_paths():
    df = pd.read_csv(Path(__file__).parent / "structure_tests.csv")
    for i, row in df.iterrows():
        facets = dict(**row)
        del facets["path"]
        yield facets, row.path


@pytest.mark.parametrize("facets,expected", get_testing_paths())
def test_build_path_from_schema(facets, expected):
    root = Path("/")
    out = build_path_from_schema(facets, root, category="raw")

    assert out is not None, "Path not build, invalid facets"
    assert str(out.relative_to(root)) == expected
