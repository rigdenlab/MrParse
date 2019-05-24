from enum import Enum
import json
import logging
import os

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
            prefix = '\033[1;{}m'.format(LogColors[record.levelname].value)
            postfix = '\033[{}m'.format(LogColors["DEFAULT"].value)
            record.msg = os.linesep.join([prefix + msg + postfix for msg in str(record.msg).splitlines()])
        return logging.Formatter.format(self, record)
 
def setup_logging():
    THIS_DIR = os.path.abspath(os.path.dirname(__file__))
    logging_json = os.path.join(THIS_DIR, 'logging.json')
    with open(logging_json, 'rt') as f:
        config = json.load(f)
    logging.config.dictConfig(config)
    return logging.getLogger()
