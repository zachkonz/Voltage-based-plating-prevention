#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Started on Thu May 11 10:06:10 2023

Author: Zachary M. Konz
This script generates Manuscript Figs. 1 and S1. First, experiment-model comparisons
are plotted for irrev. Li plating accumulation at 3 unique rates, using data from
Konz et al. 'High-throughput Li plating...' Nat. Energy (2023). Then, model results
are displayed for varying the initial charging temperature to illustrate reasonable
li plating temperature response. 
"""
import matplotlib.pyplot as plt
import pandas as pd
import utils

#%% Fig. 1: experimental vs. modeled irrev. Li data

#Experimental specs. from Konz et al., Nature Energy
n_cells = 12
fullCap_expt = 4.3 #mAh, meas., avg
grCap_expt = 5.23 #mAh, meas., avg
grCap_theoretical = 3.35 #mAh/cm2
fullCap_theoretical = 2.8 #mAh/cm2
grInitial_x = 0.015 #inital lithiation

#Each cell contains a set of SOC, irrev. Li plating columns. All cell SOC columns preceed irrev. Li 
#Note SOC here is defined as percent of full-cell capacity, fullCap_expt
df_Expt = pd.read_csv('ExperimentalFullCellsIrrevLiData_forModelComparison.csv')
colNames = df_Expt.columns.tolist()
SOC_cols = list([col for col in colNames if 'SOC' in col])
Li_cols = list([col for col in colNames if 'Li' in col])

#separate into 2 DataFrames, one for SOC, the other for irrevLi, to make indexing easier
df_SOC = df_Expt[SOC_cols]
df_irrevLi = df_Expt[Li_cols]

#Plot settings
color375C = [0, 0.38, 0.07]
color5C = [0.4039, 0.098, 0.1098]
color75C = [0.1686, 0, 0.7216]
fig, ax = plt.subplots(nrows=1,ncols=1,sharex=True)
fig.set_size_inches(4, 3.5)
plt.rcParams.update({'font.size': 12})

for i in range(n_cells):
    label = Li_cols[i]
    
    #set color and marker style according to C-rate
    if '3.75C' in label:
        color = color375C; m = 's';
    elif ' 5C' in label:
        color = color5C; m = '^';
    elif ' 7.5C' in label:
        color = color75C; m = 'o';
    
    indPlot = df_SOC.iloc[:,i] !=0 #cells cover different SOC ranges, so exclude 0 values
    # xPlot = df_SOC.loc[indPlot,SOC_cols[i]]*fullCap_expt/grCap_expt+grInitial_x #calculate Gr lithiation from full-cell SOC
    xPlot = df_SOC.loc[indPlot,SOC_cols[i]]
    ax.plot(xPlot,df_irrevLi.loc[indPlot,Li_cols[i]],('-'+m),linewidth=1,
            c=color, markerfacecolor='w',markeredgecolor=color)
    
ax.set_ylabel('Irreversible Li (%)')
ax.set_ylim([-0.02,0.79])
ax.set_xlim([0.05,0.95])
ax.set_xlabel('SOC, full-cell')

#Add the simulated data to the plot
paramNames = ['Gr_i0','C-rate','Gr_porosity','p_anode','p_sep']
outputNames = ['t', 'T', 'V', 'revLi', 'irrevLi',
            'capacity','rate','simulation #']

fileName = 'compareExptModel_calibration_rates_porosity_i0_brugg_simulations.txt'

colNames = paramNames+ outputNames
df_alldata, df_params, df_outputs = utils.createDfFromTxt(fileName, colNames)

#Note: the experimental capacity is 2.7 mAh/cm2, whereas the model
#theoretical capacity used to define the model C-rates is 2.8 mAh/cm2
#Given the slight discrepancy, the model C-rate values are slightly lower to
#achieve the same current density. 3.6C, 4.8C, 7.2C
rate_c = {3.6:color375C, 4.8:color5C, 7.2:color75C}
rate_m = {3.6:'s', 4.8:'^', 7.2:'o'}
ls = '--'


#Plot specific conditions
paramNames = ['Gr_i0','C-rate','Gr_porosity','p_anode','p_sep']
cond_eps = df_params['Gr_porosity']==0.34
cond_i0 = df_params['Gr_i0']==0.6
cond_p_an = df_params['p_anode']==2.0
cond_p_sep = df_params['p_sep']==1.8
df_params_toPlot = df_params[cond_eps & cond_i0 & cond_p_an & cond_p_sep]

for i in df_params_toPlot.index.tolist():
        rate_i = df_params['C-rate'][i]
        i0_i = df_params['Gr_i0'][i]
        p_an_i = df_params['p_anode'][i]
      
        df_outputs_i = df_outputs[df_outputs['simulation #']==i]
        irrevLi_i = df_outputs_i['irrevLi_pct']
        SOC_i = df_outputs_i['SOC_full']
        print(max(df_outputs_i['T']))#get maximum simulated temperature rise
        ax.plot(SOC_i,irrevLi_i,linestyle=ls,linewidth=3,c=rate_c[rate_i])

#add dummy lines to create legend
#7.5C, 5C, 3.75C, Model
dummy_x = [-1,-1]; dummy_y = [-0.5,-0.5];
lw_legend = 1;
ax.plot(dummy_x,dummy_y,'-o',linewidth=lw_legend,c=color75C,markeredgecolor=color75C,markerfacecolor='w',label='7.5C')
ax.plot(dummy_x,dummy_y,'-^',linewidth=lw_legend,c=color5C,markeredgecolor=color5C,markerfacecolor='w',label='5C')
ax.plot(dummy_x,dummy_y,'-s',linewidth=lw_legend,c=color375C,markeredgecolor=color375C,markerfacecolor='w',label='3.75C')
ax.plot(dummy_x,dummy_y,linestyle=ls,linewidth=3,c='k',label='model')
ax.legend(edgecolor='w',fontsize=11)

# fileName = 'Fig1_v1_exptModelCompare_24May2023.svg'
# fig.savefig(fileName,bbox_inches = 'tight', dpi= 300)
 
#%% Fig. S1: Li plating response to temperature variation
paramNames = ['Temp (K)','C-rate']
outputNames = ['t', 'T', 'V', 'revLi', 'irrevLi',
            'capacity','rate','simulation #']

fileName = 'modelTemperatureDemonstration'
colNames = paramNames+ outputNames
df_alldata, df_params, df_outputs = utils.createDfFromTxt(fileName, colNames)

fig, ax = plt.subplots(nrows=1,ncols=1,sharex=True)
fig.set_size_inches(4, 3.5)
plt.rcParams.update({'font.size': 12})

for i in range(len(df_params)):
    df_outputs_i = df_outputs[df_outputs['simulation #']==i]
    irrevLi_i = df_outputs_i['irrevLi_pct']
    SOC_i = df_outputs_i['SOC_full']*fullCap_theoretical/grCap_theoretical+grInitial_x
    # SOC_i = df_outputs_i['SOC_full']
    ax.plot(SOC_i,irrevLi_i,linewidth=2)

ax.set_ylabel('Irreversible Li (%)')
ax.set_ylim([-0.02,0.62])
ax.set_xlim([0.01,0.95])

# ax.set_xlabel('SOC, full-cell')
ax.set_xlabel('x in Li$_{x}$C$_{6}$')
legend_list = ['20 C','30 C','40 C','50 C','60 C']
ax.legend(legend_list,edgecolor='w',fontsize=10.5)
# fileName = 'FigS1_v2_modelTemperatureResponse_GrLithiation_24May2023.svg'
# fig.savefig(fileName,bbox_inches = 'tight', dpi= 300)
