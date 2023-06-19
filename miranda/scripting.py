"""Scripting Helpers module."""
from __future__ import annotations

import pathlib
import sys
from datetime import datetime as dt

_CONSOLE_FORMAT = "%(asctime)s [%(levelname)s] %(name)s: %(message)s"
_LOGFILE_FORMAT = "%(asctime)s: [%(levelname)s]: %(filename)s(%(funcName)s:%(lineno)s) >>> %(message)s"

__all__ = ["LOGGING_CONFIG"]

LOGGING_CONFIG = {
    "version": 1,
    "disable_existing_loggers": False,
    "formatters": {
        "standard": {"format": _CONSOLE_FORMAT},
        "logfile": {"format": _LOGFILE_FORMAT},
    },
    "handlers": {
        "default": {
            "level": "INFO",
            "formatter": "standard",
            "class": "logging.StreamHandler",
            "stream": sys.stdout,
        },
        "rotated_file": {
            "class": "logging.handlers.RotatingFileHandler",
            "formatter": "logfile",
            "filename": f"{dt.now().strftime('%Y%m%d')}_{pathlib.Path(sys.argv[0]).stem}.log",
            "maxBytes": 2 * 1024 * 1024,  # 2 MB
            "backupCount": 10,
        },
        "standard_file": {
            "class": "logging.FileHandler",
            "formatter": "logfile",
            "filename": f"{dt.now().strftime('%Y%m%d')}_{pathlib.Path(sys.argv[0]).stem}.log",
        },
    },
    "loggers": {
        "": {  # root logger
            "handlers": ["default", "standard_file"],
            "level": "INFO",
            "propagate": False,
        },
        "miranda": {
            "handlers": ["default", "rotated_file"],
            "level": "DEBUG",
            "propagate": False,
        },
        "__main__": {  # if __name__ == '__main__'
            "handlers": ["default"],
            "level": "INFO",
            "propagate": False,
        },
    },
}
