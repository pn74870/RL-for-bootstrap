# -*- coding: utf-8 -*-
"""
Created on Mon Jul  8 16:10:41 2024

@author: User
"""

import os, json
from Where_execute import Where_execute

root_dir = Where_execute()
def get_allowed_2states():
    file_name = 'LP_data/allowed_2delta.json'
    with open(os.path.join(root_dir, file_name), 'r') as file:
        combined_data = json.load(file)
    
    return [eval(key)['delta'] for key in combined_data.keys()]

