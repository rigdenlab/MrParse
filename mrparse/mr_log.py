from enum import Enum
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
 
