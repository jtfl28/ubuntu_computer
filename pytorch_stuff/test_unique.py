#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 29 13:46:06 2019

@author: gemsec-user
"""

import pandas as pd

#table = pd.read_csv('allpeps.txt')
with open('allpeps.txt', 'r') as g:
    table = g.readline().split(',')
    
#table = pd.DataFrame(table)
new_table = []
bad = ['X', 'Z', 'B', 'O']
for i in table:
    is_bad = False
    for j in bad:
        if j in i:
            is_bad = True
    if not is_bad:
        new_table.append(i)
table = new_table   
    
data = []
num = []
for i in range(len(table[0])):
    print(i)
    data.append({})
    for j in range(len(table)):
        if table[j][0:i+1] not in data[i]:
            data[i][(table[j][0:i+1])] = 1
        else:
            data[i][(table[j][0:i+1])] += 1
#        if len(data[i].keys()) == 20**(i+1):
#            break
    num.append((len(data[i].keys())/(20**(i+1)), len(data[i].keys())))
        