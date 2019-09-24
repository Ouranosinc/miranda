from datetime import datetime as dt

from miranda.utils import MiB

_CONSOLE_FORMAT = "%(asctime)s [%(levelname)s] %(name)s: %(message)s"
_LOGFILE_FORMAT = "%(asctime)s: [%(levelname)s]: %(filename)s(%(funcName)s:%(lineno)s) >>> %(message)s"


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
            "stream": "ext://sys.stdout",  # Default is stderr
        },
        "rotated_file": {
            "class": "logging.handlers.RotatingFileHandler",
            "formatter": "logfile",
            "filename": "{}_miranda.log".format(dt.now().strftime("%Y%m%d")),
            "maxBytes": 2 * MiB,
            "backupCount": 10,
        },
        "standard_file": {
            "class": "logging.FileHandler",
            "formatter": "logfile",
            "filename": "{}_miranda.log".format(dt.now().strftime("%Y%m%d")),
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
            "level": "WARNING",
            "propagate": False,
        },
    },
}
