# -*- coding: utf-8 -*-
"""
Created on Thu Jun 27 17:26:55 2019

@author: jtfl2
"""

import cleaning
import os
d = os.getcwd()

for i in os.listdir(d + '/pdb'):
    cleaning.cleanATOM(d + '/pdb/' + i, out_file= d + '/cleaned/' + i, ext = '.pdb')