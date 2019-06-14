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


def make_workdir(dir_name_stem='mrparse'):
    """Make a work directory rooted at run_dir and return its path

    Parameters
    ----------
    dir_name_stem : str
       name to use as stem of working directory name

    Returns
    -------
    work_dir : str
       The path to the working directory

    """
    MAX_WORKDIRS = 100
    run_dir = os.getcwd()
    run_inc = 0
    while True:
        dname = "{}_{}".format(dir_name_stem, run_inc)
        work_dir = os.path.join(run_dir, dname)
        if not os.path.exists(work_dir):
            break
        run_inc += 1
        if run_inc > MAX_WORKDIRS:
            raise RuntimeError("Too many work directories! {0}".format(work_dir))
    if os.path.exists(work_dir):
        raise RuntimeError("There is an existing work directory: {0}\n"
                           "Please delete/move it aside.".format(work_dir))
    os.mkdir(work_dir)
    return work_dir


def now():
    return datetime.datetime.now().strftime("%d/%m/%y %H:%M:%S")


def run_cmd(cmd, env={'PYTHONPATH':''}):
    """Should replace with pyjob
    
    Always set PYTHONPATH to null so processes don't inherit our environment
    """
    logger.debug("Running cmd: %s", " ".join(cmd))
    optd = { 'stderr': subprocess.STDOUT }
    if env:
        optd['env'] = env
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
