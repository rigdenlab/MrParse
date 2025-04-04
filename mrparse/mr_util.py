"""
Created on 14 Dec 2018

@author: jmht
"""
import copy
import datetime
import logging
import os
from pathlib import Path
import subprocess
import sys

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
    return fpath and Path(fpath).exists() and os.access(fpath, os.X_OK)


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
    max_work_dirs = 100
    run_dir = Path.cwd()
    run_inc = 0
    while True:
        dname = f"{dir_name_stem}_{run_inc}"
        work_dir = run_dir.joinpath(dname)
        if not work_dir.exists():
            break
        run_inc += 1
        if run_inc > max_work_dirs:
            raise RuntimeError(f"Too many work directories! {work_dir}")
    if work_dir.exists():
        raise RuntimeError(f"There is an existing work directory: {work_dir}\n"
                           f"Please delete/move it aside.")
    work_dir.mkdir()
    return str(work_dir)


def now():
    return datetime.datetime.now().strftime("%d/%m/%y %H:%M:%S")


def run_cmd(cmd):
    """Should replace with pyjob
    
    Always set PYTHONPATH to null so processes don't inherit our environment
    This needs some thinking about
    """
    logger.debug("Running cmd: %s", " ".join(cmd))
    optd = {'stderr': subprocess.STDOUT}
    pythonpath = 'PYTHONPATH'
    if pythonpath in os.environ:
        env = copy.copy(os.environ)
        env.pop(pythonpath)
        optd['env'] = env
    optd['encoding'] = 'utf-8'
    try:
        out = subprocess.check_output(cmd, **optd)
    except Exception as e:
        logger.debug("Error submitting cmd %s: %s", cmd, e)
        logger.debug("Traceback is:", exc_info=sys.exc_info())
        logger.debug("Output from job is: %s", e)
        raise (e)
    logger.debug("%s got output: %s", cmd, out)
    return out
