#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 26 11:17:24 2023

@author: Zkonz
"""

import pandas as pd

#%% #First, use an example formatted Txt file exported by COMSOL to determine the units of each parameter
fileName = 'exampleTxtFileCOMSOLformatParametricSweep.txt'
f = open(fileName)
lines=f.readlines()
paramNames = []
units = []
values = []
for x in lines:
    paramNames.append(x.split(' ')[0])
    values.append(x.split(' ')[1])
    units.append(x.split(' ')[2].split('\n')[0])

#%% From CSV file generated with generateTempCurrentStepParams.py, create Txt file in COMSOL format
df = pd.read_csv('params_noAging_1000.csv')

#drop first column, the dataframe duplicated index
df = df.drop(df.columns[0],axis=1)

#For each parameter, create a list of all values separated by commas
paramNames = df.columns.tolist()
n = len(df)
values = []
for i in range(len(paramNames)):
    list_values_i = df[paramNames[i]].tolist()
    values_i = ''
    
    for j in range(len(df)):
        values_i = values_i + str(list_values_i[j])+','
    
    values.append(values_i[:-1])
    
#%% Write parameters file with formatting observed in example file
fileName = 'exampleCOMSOLFormatted_fromCSV.txt'
with open((fileName[:-4]+'.txt'),'w') as f:
    for i in range(len(paramNames)):
        line_i = paramNames[i] + ' ' + values[i] + ' ' + units[i]
        f.write(line_i)
        f.write('\n') #go to next line
