from __future__ import annotations

import numpy as np


class TestRDRSConversion:
    def test_get_drop_vars(self, multivariable_dataset, tmp_path):
        from miranda.convert.eccc_rdrs import _get_drop_vars

        ds = multivariable_dataset(
            values=np.linspace(0, 1, 10),
            variables=["tmax", "tmin", "tdps", "evspsblpot", "pr"],
        )
        file = tmp_path / "test_rdrs.nc"
        ds.to_netcdf(file)

        keep_vars = ["tmax"]

        dropped = _get_drop_vars(file, keep_vars=keep_vars)

        assert set(dropped) == {
            "evspsblpot",
            "tmin",
            "pr",
            "tdps",
        }, "The dropped variables do not match the expected ones."
