from __future__ import annotations

import numpy as np
import pytest
import xarray as xr


try:
    import regionmask
except ImportError:
    regionmask = None

try:
    import clisops
except ImportError:
    clisops = None


class TestDomains:
    @pytest.mark.skipif(regionmask is None, reason="regionmask not installed")
    def test_ar6_regions(self, cassini):
        from miranda.gis._domains import add_ar6_regions

        ds = xr.open_dataset(cassini.fetch("CMIP6/snw_day_CanESM5_historical_r1i1p1f1_gn_19910101-20101231.nc"))
        ar6 = add_ar6_regions(ds)

        assert "region" in ar6.coords

    @pytest.mark.skipif(regionmask is not None, reason="regionmask is installed")
    def test_ar6_regions_no_regionmask(self, cassini):
        from miranda.gis._domains import add_ar6_regions

        ds = xr.open_dataset(cassini.fetch("CMIP6/snw_day_CanESM5_historical_r1i1p1f1_gn_19910101-20101231.nc"))

        with pytest.raises(
            ImportError,
            match="`add_ar6_regions` requires installation of the miranda GIS libraries.",
        ):
            add_ar6_regions(ds)

    @pytest.mark.parametrize(
        "domain, expected",
        [
            ("global", [90.0, -180.0, -90.0, 180.0]),
            ("AMNO", [90.0, -179.9, 10.0, -10.0]),
            ("CAN", [83.5, -141.0, 41.5, -52.5]),
            ("QC", [63.0, -80.0, 44.5, -57.0]),
            ("mtl", [45.75, -74.05, 45.3, -73.4]),
            ("fiji", "NotImplemented"),
        ],
    )
    def test_subsetting_domains(self, domain, expected):
        from miranda.gis._domains import subsetting_domains

        if expected == "NotImplemented":
            with pytest.raises(NotImplementedError):
                subsetting_domains(domain)
        else:
            result = subsetting_domains(domain)
            assert result == expected

    @pytest.mark.skipif(clisops is None, reason="clisops not installed")
    @pytest.mark.parametrize(
        "dataset",
        ["ECMWF/era5_t2m-d2m_2024.nc", "NASA/daymet_v4_prcp_annttl_hi_1987.nc"],
        ids=["era5", "daymet"],
    )
    def test_subset_domain(self, cassini, dataset):
        from miranda.gis import subset_domain

        ds = xr.open_dataset(cassini.fetch(dataset))

        if "longitude" in ds.coords and "latitude" in ds.coords:
            # Rename coordinates to lon/lat if they are named differently
            ds = ds.rename({"longitude": "lon", "latitude": "lat"})

        if dataset == "NASA/daymet_v4_prcp_annttl_hi_1987.nc":
            with pytest.raises(ValueError):
                subset_domain(ds, "mtl")
        else:
            subset = subset_domain(ds, "mtl")

            np.testing.assert_array_equal(subset.lon.values, [-74.0, -73.75, -73.5])
            np.testing.assert_array_equal(subset.lat.values, [45.75, 45.5])
