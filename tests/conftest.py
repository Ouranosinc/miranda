from __future__ import annotations

from pathlib import Path

data_path = Path.cwd()

test_data = dict()
test_data["current"] = data_path / "cmip5" / "tmax.cur.nc"
test_data["future"] = data_path / "cmip5" / "tmax.fut.nc"
test_data["observed"] = data_path / "cmip5" / "tmax.obs.nc"
