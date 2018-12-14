'''
Created on 14 Dec 2018

@author: jmht
'''

import datetime

def now():
    return datetime.datetime.now().strftime("%d/%m/%y %H:%M:%S")
