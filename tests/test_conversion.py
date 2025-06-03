from __future__ import annotations

import numpy as np

from miranda.convert import aggregations_possible


def test_possible_aggregations(multivariable_dataset):
    """Test the aggregations_possible function."""
    variables = [
        "tasmax",
        "tasmin",
        "tdps",
        "evspsblpot",
        "pr",
        "rdrs",
        "hurs",
        "ua",
        "va",
    ]

    values = np.random.rand(10)
    ds = multivariable_dataset(values, variables=variables)

    # Check possible aggregations for daily frequency
    result = aggregations_possible(ds, freq="day")
    assert isinstance(result, dict)
    # TODO: These might not be the final expected results, adjust as needed
    assert result == {
        "_tas": {
            "max",
            "mean",
            "min",
        },  # tas has an underline to denote it will be calculated from tasmax and tasmin
        "tasmax": {"max", "mean", "min"},
        "tasmin": {"max", "mean", "min"},
        "tdps": {"max", "mean", "min"},
        "evspsblpot": {"mean"},
        "pr": {"mean"},
        "hurs": {"max", "mean", "min"},
        "ua": {"mean"},
        "va": {"mean"},
    }
