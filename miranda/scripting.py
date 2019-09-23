import logging
from datetime import datetime as dt
from pathlib import Path
from typing import Union

_LOGFILE_FORMAT = "%(asctime)s %(name)s: %(levelname)s: %(message)s"
_CONSOLE_FORMAT = "%(filename)s: {}".format(_LOGFILE_FORMAT)


def logger_setup(
    file: Union[str, Path] = None,
    name: str = None,
    log_level: int = logging.DEBUG,
    console_level: int = logging.INFO,
) -> logging.Logger:
    if not file:
        file = "miranda"
    else:
        file = Path(file).stem
    logfile = "{}_{}.log".format(dt.now().strftime("%Y%m%d"), file)

    # Console formatting
    console_formatter = logging.Formatter(_CONSOLE_FORMAT)
    console_logger = logging.StreamHandler()
    console_logger.setLevel(console_level)
    console_logger.setFormatter(console_formatter)

    # Logfile formatting
    logfile_formatter = logging.Formatter(_LOGFILE_FORMAT)
    file_logger = logging.FileHandler(filename=logfile)
    file_logger.setLevel(log_level)
    file_logger.setFormatter(logfile_formatter)

    logger = logging.getLogger(name)
    logger.setLevel(logging.DEBUG)
    logger.addHandler(console_logger)
    logger.addHandler(file_logger)
    return logger
