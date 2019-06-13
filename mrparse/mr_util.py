"""
Created on 14 Dec 2018

@author: jmht
"""

import datetime
import logging
import os
import subprocess
import sys

from pyjob.script import EXE_EXT


PYTHONVERSION = sys.version_info[0]


logger = logging.getLogger(__name__)


def is_exe(fpath):
    """Check if an executable exists

    Parameters
    ----------
    fpath : str
       The path to the executable

    Returns
    -------
    bool

    """
    return fpath and os.path.exists(fpath) and os.access(fpath, os.X_OK)


def now():
    return datetime.datetime.now().strftime("%d/%m/%y %H:%M:%S")


def run_cmd(cmd):
    """Should replace with pyjob"""
    logger.debug("Running cmd: %s", " ".join(cmd))
    optd = { 'stderr': subprocess.STDOUT }
    if PYTHONVERSION > 2:
        optd['encoding'] = 'utf-8'
    try:
        out = subprocess.check_output(cmd, **optd)
    except Exception as e:
        logger.debug("Error submitting cmd %s: %s", cmd, e)
        logger.debug("Traceback is:", exc_info=sys.exc_info())
        logger.debug("Output from job is: %s", e.output)
        raise(e)
    logger.debug("%s got output: %s", cmd, out)
    return out
