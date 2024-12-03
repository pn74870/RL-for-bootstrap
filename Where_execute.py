# -*- coding: utf-8 -*-
"""
Created on Fri Aug 25 16:01:18 2023

@author: User
"""
import os
def Where_execute():
    path = '/Users/justinas/Documents/python/lp ml/'
    return path
    if os.path.isdir(path):
        return path
    else:
        return '/home/string-3/lp-env/'
