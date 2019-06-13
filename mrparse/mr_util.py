"""
Created on 14 Dec 2018

@author: jmht
"""

import datetime
import os


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
