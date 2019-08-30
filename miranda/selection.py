import logging
from datetime import date
from datetime import datetime as dt
from pathlib import Path
from typing import Union

import fabric

from miranda.utils import creation_date


def move_by_date(
    source: Union[Path, str],
    target: Union[Path, str],
    server: Union[Path, str],
    username: str,
    password: str,
    year: int,
    month: int,
    day: int,
    datetime: date = None,
) -> None:

    date_selected = date(year, month, day) or datetime

    nc_files = Path(source).glob("**/*.nc")
    nc_files = list(nc_files)
    nc_files.sort()

    connection = fabric.Connection(
        host=server, user=username, connect_kwargs=dict(password=password)
    )

    with connection as context:
        for file in nc_files:
            if creation_date(file) == date_selected:
                logging.info(
                    "{}: Moving file {} to {}".format(
                        dt.now().strftime("%Y-%m-%d %X"), file, target
                    )
                )

                if not Path(target, file.stem).exists():
                    context.get(file, target)
    return
