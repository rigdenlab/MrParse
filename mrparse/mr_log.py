from enum import Enum
import json
import logging
import logging.config
import os
from pathlib import Path


class LogColors(Enum):
    """Color container for log messages"""
    CRITICAL = 31
    DEBUG = 34
    DEFAULT = 0
    ERROR = 31
    WARNING = 33


class LogColorFormatter(logging.Formatter):
    """Formatter for log messages"""
    def format(self, record):   
        if record.levelname in LogColors.__members__:
            prefix = f'\033[1;{LogColors[record.levelname].value}m'
            postfix = f'\033[{LogColors["DEFAULT"].value}m'
            record.msg = os.linesep.join([prefix + msg + postfix for msg in str(record.msg).splitlines()])
        return logging.Formatter.format(self, record)
 
def setup_logging():
    THIS_DIR = Path(__file__).resolve().parent
    logging_json = Path(THIS_DIR, 'logging.json')
    with open(logging_json, 'rt') as f:
        config = json.load(f)
    logging.config.dictConfig(config)
    return logging.getLogger()
