# -*- coding: utf-8 -*-
"""
This version started on Mon Apr 10 2023

@author: Zachary M Konz

Analyze multistep fast-charge simulation results of Li plating for varied
temperatures, currents, SOC, and extent of aging. Code used to produce Manuscript
Figures 4-6, S4-S7.

The simulations have specified temperature, current profiles and aging parameters
according to generateTempCurrentStepParams.py

This script utilizes 4 data files that are assumed to be in the same project folder:
    1. "fullOCVexptCo10extended.csv" - contains experimental OCV vs. SOC data similar to model
    2. "PlatingBoundaryCurves.csv" - contains V vs. SOC data columns for each boundary curve
        -> Note, these were created by sketching the curves along the data points 
        in MS Powerpoint and then using webplotdigitizer
    3. "multistepProtocol_noAging_1000simulations.txt" - contains COMSOL output
    4. "multistepProtocol_correlatedAgingParamsByHealth_2000simulations.txt" - contains COMSOl output

"""
#%% Import data and packages
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import seaborn as sns
import warnings
import utils
from scipy.interpolate import splrep, BSpline

warnings.filterwarnings("ignore")

OCV_data = np.loadtxt("fullOCVexptCo10extended.csv", delimiter = ',')# C/10 Gr|NMC532 charging curve
OCV_x = OCV_data[:,0].ravel();
OCV_V = OCV_data[:,1].ravel();

#Upload Voltage Boundary curves for SOH 0.85,0.9,0.95,0.98,1.0 
#data file contains 10 columns [SOC_0.85, V_0.85, .... , SOC_1, V_1]
df_boundary_data = pd.read_csv("PlatingBoundaryCurves.csv")

paramNames = ['t1','t2','t3','t4','t5','t6','t7','t8',
              'rate1','rate2','rate3','rate4',
              'ramp1','ramp2','ramp3','ramp4','ramp5','ramp6','ramp7','ramp8',
              'T0','T1','T2','T3','T4','T5','T6','T7',
             'T target', 'iSOC']

agingParamNames = ['an_expansion_um','ca_expansion_um','an_void','ca_void','an_i0_factor',
            'ca_i0_factor','cond_factor','an_LiLoss','an_LAM','ca_LAM','SOH']

outputNames = ['t', 'T', 'V', 'revLi', 'irrevLi',
            'capacity','rate','simulation #']

#%% Import simulations v1, no aging
fileName = "multistepProtocol_noAging_1000simulations.txt"
colNames = paramNames+ agingParamNames + outputNames
df_alldata_BOL, df_params_BOL, df_outputs_BOL = utils.createDfFromTxt(fileName, colNames)

#%% Import simulations v2, SOH-correlated aging
fileName ="multistepProtocol_correlatedAgingParamsByHealth_2000simulations.txt"
colNames = paramNames + agingParamNames + outputNames
df_alldata_SOH, df_params_SOH, df_outputs_SOH = utils.createDfFromTxt(fileName, colNames)

#%% Calculate plating onsets using a threshold, store data for later
onset_threshold = 0.01 #irreversible Li as percentage of graphite capacity, %
df_outputs = df_outputs_SOH
df_params = df_params_SOH

df_onsets, V_all, SOC_all = utils.parseData_calcOnsets(df_outputs, onset_threshold)
df_simInfo_original = pd.concat([df_params,df_onsets],axis=1)
df_simInfo_original['iSOC_err'] = 0 #initialize SOC_err to 0
#%% Functions to filter and preprocess data
def filterByOnset(df_simInfo, cap_min, SOC_max, V_min, V_max, sameRate=False):
    """
        This function returns the indices of df_simInfo that satisfy the provided
        criteria for Li plating onset conditions
    """ 
    cond1 = df_simInfo['Onset Cap']>cap_min
    cond2 = df_simInfo['Onset SOC']<SOC_max
    cond3 = df_simInfo['Onset Voltage']>V_min
    cond4 = df_simInfo['Onset Voltage']<V_max
    cond5 = True
    
    if sameRate:
        cond5 = df_simInfo['Onset rate']==df_simInfo['rate_4']
        #add condition for no rate change within 
    
    filtered_df_onsets = df_simInfo[cond1 & cond2 & cond3 & cond4 & cond5]
    
    return filtered_df_onsets.index

#%% Fig. 6b-e. Visualize data plating onsets, colored by SOH, with V(SOC) boundaries
df_simInfo = df_simInfo_original.copy()

#first filter plating onsets - primarily to remove simulations that do not induce Li plating
ind_filtered = filterByOnset(df_simInfo, cap_min=0.0, SOC_max=0.95, V_min=3.8, V_max=4.39)

#Plotting options and Settings
plotBySOH = True #if set to True, then plot the onset stars
plotBoundary = True #plots voltage boundary curve
onSamePlot = False #plots all data or boundary curves on same plot, else individual plots by SOH
plotOCV = False #plots the experimental OCV curve for comparison w/ V(SOC) boundary
saveFigs = False #saves figures as SVG if True
SOH_toPlot = [0.85,0.90,0.95,0.98]#specify which SOH to plot, of 0.85,0.9,0.95,0.98
SOC_max_confident = 0.95 #upper x limit of boundary curve plotted
SOC_min_confident = 0.1 #lower x limit of boundary curve plotted
plt.rcParams.update({'font.size': 13})
fs=13
lw = 2

#Plotting
if onSamePlot:
    #create figure
    fig, ax = plt.subplots(nrows=1,ncols=1)
    fig.set_size_inches(6, 3.5)
    ax.set_ylabel('Voltage (V)', fontsize=fs)
    ax.set_yticks([3.8,4,4.2,4.4])
    ax.set_xlabel('State of Charge', fontsize=fs)
    ax.set_ylim([3.78,4.41])
    ax.set_xlim([-0.02,1])

for SOH_i in SOH_toPlot:

    if not onSamePlot: #create individual figure for each SOH
        fig, ax = plt.subplots(nrows=1,ncols=1)
        fig.set_size_inches(6, 3)
        ax.set_ylabel('Voltage (V)', fontsize=fs)
        ax.set_xlabel('State of Charge', fontsize=fs)
        ax.set_ylim([4.03,4.41])
        ax.set_xlim([-0.02,1])
    
    if plotBySOH:
        df_filtered = df_simInfo.loc[ind_filtered,:]
        data_filtered_SOH = df_filtered[df_filtered['SOH']==SOH_i]
        # print(len(data_filtered_SOH))
        sns.scatterplot(data = data_filtered_SOH, x = 'Onset SOC', y = 'Onset Voltage',
                        ax=ax, hue = 'SOH', hue_norm=(0.82,1.02), palette = 'coolwarm_r', legend=False,
                        size = 'Onset rate', sizes=(100,150),size_norm=(0,8), marker='*',
                        edgecolor='k', linewidth=0.5)

    if plotOCV:
        ax.plot(OCV_x/max(OCV_x),OCV_V,linestyle=':',linewidth=lw,c='k',label='OCV')

    if plotBoundary:
        #Plot boundaries assigned from WebplotDigitizer:
        
        #calculate boundary curve
        SOC_label='SOC_' + str(SOH_i)
        V_label = 'V_' + str(SOH_i)
    
        df_boundary_nonzero = df_boundary_data[df_boundary_data[V_label]>0]
        SOC_b = df_boundary_nonzero[SOC_label]
        V_b = df_boundary_nonzero[V_label]
        
        #need to calculate spline
        minV_SOC_spline = splrep(SOC_b,V_b, s=0.0005)
        
        cmap = matplotlib.cm.get_cmap('coolwarm_r')
        norm = matplotlib.colors.Normalize(vmin=0.82,vmax=1.02)
        c = cmap(norm(SOH_i))
        
        #plot boundary
        SOC_test = np.linspace(SOC_min_confident,SOC_max_confident,50)
        ax.plot(SOC_test, BSpline(*minV_SOC_spline)(SOC_test), color=c,
                linewidth=lw+0.5,label=None)
    if saveFigs: 
        fileName = 'onsetStars_SOH_all_wBoundary'+ str(SOH_i)+'.svg'
        fig.savefig(fileName,bbox_inches = 'tight', dpi= 300)
        
#%% Functions to calculate voltage curve intersection with boundary
def addVandCap_boundaryIntersect(df_simInfo, df_outputs, boundary_spline, y_offset=0,x_offset=0):
    """
        This function finds the first intersection (lowest capacity) of the simulated
        V(SOC) voltage profile with the voltage boundary spline. These are termed
        Voltage and Capacity 'cutoffs', and are added for each simulation to the
        df_simInfo DataFrame
    """ 
    
    list_Vcutoffs_all = []
    list_Capcutoffs_all = []
    colNames = ['V_cutoff','Cap_cutoff']
    
    for i in df_simInfo.index:
        df_outputs_i = df_outputs[df_outputs['simulation #']==i]
        V_i = np.array(df_outputs_i['V'])
        SOC_i = np.array(df_outputs_i['chargeCap'] + df_simInfo['iSOC'][i] + df_simInfo['iSOC_err'][i])+x_offset
        V_threshold = BSpline(*boundary_spline)(SOC_i)-y_offset
        SOC_cutoff_i, V_cutoff_i = getFirstIntersection(SOC_i,V_i,V_threshold)
        list_Vcutoffs_all.append(V_cutoff_i)
        list_Capcutoffs_all.append(SOC_cutoff_i-df_simInfo['iSOC'][i])
    
    #add columns to df_simInfo
    data = np.array([list_Vcutoffs_all,list_Capcutoffs_all]).T
    
    for i in range(len(colNames)):
        df_simInfo[colNames[i]] = data[:,i]

    return df_simInfo

def getFirstIntersection(SOC,V_meas,V_threshold):
    """
        This function is called by addVandCap_boundaryIntersect(). The inputs
        are an SOC vector and associated model calculated voltage (V_meas) and 
        voltage boundary values for that SOC
    """ 
    #find first point where V_meas exceeds V_threshold
    ind = (V_meas>V_threshold).nonzero()[0][0]
    
    #use that and previous point to calculate point of intersection with threshold curve
    V_meas2 = V_meas[ind]
    V_meas1 = V_meas[ind-1]
    V_spl2 = V_threshold[ind]
    V_spl1 = V_threshold[ind-1]
    SOC_1 = SOC[ind-1]
    SOC_2 = SOC[ind]
    
    SOC, V = calcIntersect(SOC_1,SOC_2,V_meas1,V_meas2,V_spl1,V_spl2)
        
    return SOC, V

def calcIntersect(x1,x2,y1i,y2i,y1ii,y2ii):
    """
        This function calculates the intersection of 2 line segments, i and ii,
        that share x-coordinates 1 and 2. It is called by getFirstIntersection()
    """ 
    
    m_i = (y2i-y1i)/(x2-x1)
    b_i = y1i-m_i*x1
    
    m_ii = (y2ii-y1ii)/(x2-x1)
    b_ii = y1ii-m_ii*x1
    
    x = (b_i-b_ii)/(m_ii-m_i)
    y = m_i*x+b_i
    return x, y

#%% For beginning-of-life (no aging) data, calculate the difference between 
#boundary intersect and plating onset, also calculate % charge completion

#first, get BOL data and calculate plating onsets
df_outputs = df_outputs_BOL
df_params = df_params_BOL
df_onsets, V_all, SOC_all = utils.parseData_calcOnsets(df_outputs, onset_threshold)
df_simInfo_original = pd.concat([df_params,df_onsets],axis=1)
df_simInfo_original['iSOC_err'] = 0 #initialize SOC_err to 0

df_simInfo = df_simInfo_original.copy()

#need to calculate boundary spline
SOC_label='SOC_1' #1 is for full health/SOH, no aging
V_label = 'V_1'
df_boundary_nonzero = df_boundary_data[df_boundary_data[V_label]>0]
SOC_b = df_boundary_nonzero[SOC_label]
V_b = df_boundary_nonzero[V_label]
boundary_spline_BOL = splrep(SOC_b,V_b, s=0.0005)

#include only simulations with plating
ind_i = filterByOnset(df_simInfo, cap_min=0.0, SOC_max=0.95, V_min=3.8, V_max=4.39)
df_simInfo_plating = df_simInfo.loc[ind_i,:]

#add voltage and capacity intersections 'cutoffs'
df_simInfo_plating = addVandCap_boundaryIntersect(df_simInfo_plating, df_outputs,boundary_spline_BOL)

#calculate SOC to onset after boundary and Pct charge to boundary
df_simInfo_plating['SOC to onset after boundary'] = df_simInfo_plating['Onset Cap']-df_simInfo_plating['Cap_cutoff']
df_simInfo_plating['Pct charge to boundary'] = 100*df_simInfo_plating['Cap_cutoff']/df_simInfo_plating['Onset Cap']
df_simInfo_plating['SOH']=1

#%% Fig. 5g-h. Plot voltage boundary effectiveness histograms for no aging
plotDeltaSOC = False #Toggle this. If False, then it will plot the % charge acceptance

plt.rcParams.update({'font.size': 30})
sns.set_style('white')

if plotDeltaSOC:
    x_plot = 'SOC to onset after boundary'
    bw = 0.007
    
    sns.displot(df_simInfo_plating, x=x_plot, hue = 'SOH',
                hue_norm=(0.82,1.02), palette = 'coolwarm_r',alpha=1, legend=False,
                binwidth=bw, height=4, aspect=2)
    
    plt.xlim([-0.01,0.21])
    plt.xlabel('$\Delta$SOC to onset after boundary')
    
else:
    x_plot = 'Pct charge to boundary'
    bw = 2.5
    
    sns.displot(df_simInfo_plating, x=x_plot, hue = 'SOH',
                hue_norm=(0.82,1.02), palette = 'coolwarm_r',alpha=1, legend=False,
                binwidth=bw, height=4, aspect=2)

    plt.xlim([-1,101])
    plt.xticks([0,50,100])
    plt.xlabel('% Charge completion')


print(df_simInfo_plating[['SOC to onset after boundary']].describe())
print(df_simInfo_plating[['Pct charge to boundary']].describe())

#%% For aging data, plot the distribution of onset-boundary deltaSOC by health index
#first, get aging data
df_outputs = df_outputs_SOH
df_params = df_params_SOH

df_onsets, V_all, SOC_all = utils.parseData_calcOnsets(df_outputs, onset_threshold)
df_simInfo_original = pd.concat([df_params,df_onsets],axis=1)
df_simInfo_original['iSOC_err'] = 0 #initialize SOC_err to 0

#get data
df_simInfo = df_simInfo_original.copy()
SOH_suffix = ['0.98','0.95','0.9','0.85']

df_list_bySOH = []

for SOH_i in SOH_suffix:
    #get associated data
    ind_SOH_i = df_simInfo['SOH']==float(SOH_i)
    df_simInfo_SOH_i = df_simInfo[ind_SOH_i]
    
    #the index of df_simInfo_SOH_i specifies which examples to compute, so don't need to reduce df_outputs
    #compute the boundary
    V_label = 'V_' + SOH_i
    SOC_label = 'SOC_' + SOH_i
    
    df_boundary_nonzero = df_boundary_data[df_boundary_data[V_label]>0]
    SOC_b = df_boundary_nonzero[SOC_label]
    V_b = df_boundary_nonzero[V_label]
    boundary_spline_BOL = splrep(SOC_b,V_b, s=0.0005)

    #only include simulations with plating 
    ind_i = filterByOnset(df_simInfo_SOH_i, cap_min=0.0, SOC_max=0.95, V_min=3.8, V_max=4.39)
    df_simInfo_plating = df_simInfo_SOH_i.loc[ind_i,:]


    #add voltage and capacity intersections 'cutoffs'
    df_simInfo_plating = addVandCap_boundaryIntersect(df_simInfo_plating, df_outputs,boundary_spline_BOL)

    df_simInfo_plating['SOC to onset after boundary'] = df_simInfo_plating['Onset Cap']-df_simInfo_plating['Cap_cutoff']
    df_simInfo_plating['V change to onset after boundary'] = df_simInfo_plating['Onset Voltage']-df_simInfo_plating['V_cutoff']
    df_simInfo_plating['Pct charge to boundary'] = 100*df_simInfo_plating['Cap_cutoff']/df_simInfo_plating['Onset Cap']

    #Uncomment the lines below to get distribution statistics for each SOH
    
    # print(df_simInfo_plating[['SOC to onset after boundary']].describe())
    # print(df_simInfo_plating[['Pct charge to boundary']].describe())
    # print(df_simInfo_plating[['V change to onset after boundary']].describe())
    
    df_list_bySOH.append(df_simInfo_plating[['SOH','SOC to onset after boundary',
                                             'V change to onset after boundary','Pct charge to boundary']])

df_allDiff = pd.concat(df_list_bySOH)

#%% Fig. 6g-h. Plot multiple health distributions
plotDeltaSOC = False #Toggle this. If False, then it will plot the % charge acceptance

plt.rcParams.update({'font.size': 30})
sns.set_style('white')

if plotDeltaSOC:
    x_plot = 'SOC to onset after boundary'
    bw = 0.007
    
    sns.displot(df_allDiff, x=x_plot, hue = 'SOH',
                hue_norm=(0.82,1.02), palette = 'coolwarm_r',alpha=1, legend=False,
                binwidth=bw, height=4, aspect=2,stat='probability', common_norm=False)
    
    plt.xlim([-0.01,0.21])
    plt.xlabel('$\Delta$SOC to onset after boundary')
    
else:
    x_plot = 'Pct charge to boundary'
    bw = 2.5
    
    sns.displot(df_allDiff, x=x_plot, hue = 'SOH',
                hue_norm=(0.82,1.02), palette = 'coolwarm_r',alpha=1, legend=False,
                binwidth=bw, height=4, aspect=2,stat='probability', common_norm=False)

    plt.xlim([-1,101])
    plt.xticks([0,50,100])
    plt.xlabel('% Charge completion')

print(df_simInfo_plating[['SOC to onset after boundary']].describe())
print(df_simInfo_plating[['Pct charge to boundary']].describe())

#%% Fig. 6a. Plot gaussian distributions from which parameter values were sampled as f(SOH)
import scipy.stats as stats
import matplotlib.cm
# import math
plt.rcParams.update({'font.size': 16})
fig, ax = plt.subplots(nrows=1,ncols=1)
fig.set_size_inches(6, 3)
fs=15
lw = 2

SOH_min = 0.8
SOH_max = 1
#for each SOH value, calculate a gaussian curve and plot
soh_vals = [0.98,0.95,0.9,0.85]

for SOH_val in soh_vals:
    sd = (1-SOH_val)/3
    mean = SOH_val
    x = np.linspace(0.7,1.05,200)
    y = stats.norm.pdf(x, mean, sd)
    
    ax.plot(x, y, color='k', linewidth=0.5)    
    
    #fill color
    cmap = matplotlib.cm.get_cmap('coolwarm_r')
    SOH_norm = (SOH_val-0.82)/(1.02-0.82)
    rgba = cmap(SOH_norm)
    ax.fill_between(x,y,color=rgba)
    
ax.plot([SOH_max, SOH_max],[-100,100],linewidth=0.5,color = 'k',linestyle=':')
ax.plot([SOH_min, SOH_min],[-100,100],linewidth=0.5,color = 'k',linestyle=':') 

ax.set_xlabel('Aging parameter value', fontsize=fs)
ax.set_ylim([-4,62])
ax.set_xlim([0.78,1.02])
ax.invert_xaxis()
ax.get_yaxis().set_visible(False)
ax.set_ylabel('Probability density', fontsize=fs)

#%% Fig. 4a-d. Plot current, voltage, temperature profiles - for 4 examples only, no filtering yet
onset_threshold = 0.01 #irreversible Li as percentage of graphite capacity, %
df_outputs = df_outputs_BOL
df_params = df_params_BOL

df_onsets, V_all, SOC_all = utils.parseData_calcOnsets(df_outputs, onset_threshold)
df_simInfo_original = pd.concat([df_params,df_onsets],axis=1)
df_simInfo_original['iSOC_err'] = 0 #initialize SOC_err to 0

df_simInfo = df_simInfo_original.copy()
ind_i = df_simInfo.index
from itertools import cycle

graphiteCap = 3.35
plt.rcParams.update({'font.size': 16})
fig, ax = plt.subplots(nrows=4,ncols=1,sharex=True,layout='constrained')
fig.set_size_inches(6, 7)
lines = ["-"]
colorcycler = cycle(['b','m','k','g'])
linecycler = cycle(lines)
lw = 1
count = 1
# for i in np.arange(0,len(ind_i),1):
for i in [64,67,87,122]:   
    df_i = df_outputs[df_outputs['simulation #']==ind_i[i]]
    df_i = df_i[df_i['SOC_full']<df_simInfo['Onset SOC'][ind_i[i]]]

    ls = next(linecycler)
    c=next(colorcycler)
    ax[3].plot(df_i['SOC_full'],df_i['V'],linewidth=lw,linestyle=ls,color=c,label='Ex. '+str(count))
    ax[2].plot(df_i['SOC_full'],df_i['irrevLi_pct']*graphiteCap*1000/100,linewidth=lw,linestyle=ls,color=c)
    ax[1].plot(df_i['SOC_full'],df_i['T'],linewidth=lw,linestyle=ls,color=c)
    ax[0].plot(df_i['SOC_full'], df_i['rate'],linewidth=lw,linestyle=ls,color=c)
    print(i)

    count = count+1
    
    ax[3].scatter(df_onsets['Onset SOC'][ind_i[i]], df_onsets['Onset Voltage'][ind_i[i]],
                  marker='*', color=[1, 1, 0], s=200, edgecolors = c,zorder=2) 
    ax[2].scatter(df_onsets['Onset SOC'][ind_i[i]], max(df_i['irrevLi_pct']*graphiteCap*1000/100),
                    marker='*', color=[1, 1, 0], s=200, edgecolors = c,zorder=2) 
    

ax[0].set_ylabel('C-rate')
ax[1].set_ylabel('T (' + chr(176) +'C)')
ax[2].set_ylabel('Li ($\mu$Ah/cm$^2$)')
ax[2].set_yticks([0.1,0.2,0.3])
ax[3].set_ylabel('Voltage')
ax[3].set_xlabel('State of Charge')

ax[0].set_ylim([1.9,8.1])
ax[0].set_yticks([2,4,6,8])
ax[1].set_ylim([8,67])
ax[3].set_ylim([3.59,4.44])
ax[3].set_yticks([3.6,4,4.4])
ax[3].set_xlim([-0.01,1.01])

#%% Fig. 5a-d. Plot current, voltage, temperature profiles - all simulations with plating, with filtering
def getValsDiscreteSOC(SOC, values, discreteSOC):
    """
    Goal: save any simulation value such as T or rate at discrete SOC, to be
    used later for computing the average values across all n=500 simulations

    Input:
      SOC (series): full-cell SOC values for the entire charge simulation
      values (series): T, rate, irrev. Li, or other values for the entire charge simulation
      discreteSOC (array): the SOC at which to interpolate values (0.02, 0.04,..0.94 for example)
    
    Output:
      discreteValues (array): interpolated values, or zeros if discreteSOC>max(SOC)
    """    

    discreteValues = np.interp(discreteSOC,SOC,values,left=0, right=0)
    
    return discreteValues

df_simInfo = df_simInfo_original.copy()
ind_i = filterByOnset(df_simInfo, cap_min=0.0, SOC_max=0.95, V_min=3.8, V_max=4.39)

graphiteCap = 3.35
fig, ax = plt.subplots(nrows=4,ncols=1,sharex=True,layout='constrained')
fig.set_size_inches(6, 7)
plt.rcParams.update({'font.size': 16})
lines = ["-"]
# linecycler = cycle(lines)
lw = 1
count = 1

discreteSOC = np.linspace(0.04,0.92,num=45) #0.04, 0.04, ... , 0.92
list_rateArrays = []
list_TArrays = []

for i in np.arange(0,len(ind_i),1):  
    df_i = df_outputs[df_outputs['simulation #']==ind_i[i]]
    df_i = df_i[df_i['SOC_full']<df_simInfo['Onset SOC'][ind_i[i]]]
    
    ax[3].plot(df_i['SOC_full'],df_i['V'],linewidth=lw,linestyle=ls,label='Ex. '+str(count))
    ax[2].plot(df_i['SOC_full'],df_i['irrevLi_pct']*graphiteCap*1000/100,linewidth=lw,linestyle=ls)
    ax[1].plot(df_i['SOC_full'],df_i['T'],linewidth=lw,linestyle=ls)
    ax[0].plot(df_i['SOC_full'], df_i['rate'],linewidth=lw,linestyle=ls)
    count = count+1

    list_rateArrays.append(getValsDiscreteSOC(df_i['SOC_full'],df_i['rate'],discreteSOC))
    list_TArrays.append(getValsDiscreteSOC(df_i['SOC_full'],df_i['T'],discreteSOC))

arr_2D_rate = np.array(list_rateArrays)
avgRates = np.mean(arr_2D_rate,axis=0,where=(arr_2D_rate>0))

arr_2D_T = np.array(list_TArrays)
avgTs = np.mean(arr_2D_T,axis=0,where=(arr_2D_T>0))

ax[0].plot(discreteSOC,avgRates,linewidth=4,linestyle=':',color='k',label='avg.')
ax[1].plot(discreteSOC,avgTs,linewidth=4,linestyle=':',color='k',label='avg.')

ax[0].set_ylabel('C-rate')
ax[0].legend(edgecolor='w',fontsize=14,loc='upper right')
ax[1].set_ylabel('T (' + chr(176) +'C)')
ax[1].legend(edgecolor='w',fontsize=14,loc='lower right')
ax[2].set_ylabel('Li ($\mu$Ah/cm$^2$)')
ax[2].set_yticks([0.1,0.2,0.3])
ax[3].set_ylabel('Voltage')
ax[3].set_xlabel('State of Charge')

ax[0].set_ylim([1.9,8.2])
ax[0].set_yticks([2,4,6,8])
ax[1].set_ylim([8,64])
ax[3].set_ylim([3.56,4.44])
ax[3].set_yticks([3.6,4,4.4])
ax[3].set_xlim([-0.01,1.01])

#%% Fig. 5e. Plot Plating onset stars and boundary on plot
df_simInfo = df_simInfo_original.copy()
ind_i = filterByOnset(df_simInfo, cap_min=0.0, SOC_max=0.95, V_min=3.8, V_max=4.39)
lw=1
plt.rcParams.update({'font.size': 15})
fig, ax = plt.subplots(nrows=1,ncols=1)
fig.set_size_inches(6, 4)
SOC_max_confident=0.93
SOC_min_confident=0.12
ax.set_ylabel('Voltage (V)')
ax.set_ylim([3.40,4.44])
# ax.set_ylim([3.9,4.2])
ax.set_xlabel('State of Charge')

#Plot 4.2V horizontal line for typical cutoff voltage
# ax.plot([0,1],[4.2,4.2],linestyle='-',linewidth=lw-0.2,c='k')

ax.scatter(df_simInfo['Onset SOC'][ind_i], df_simInfo['Onset Voltage'][ind_i], marker='*', color=[1, 1, 0], s=100, 
            edgecolors = 'k',linewidth=0.5,zorder=2 ,label = "Plating Onset")      

#Plot OCV curve
ax.plot(OCV_x/max(OCV_x)*0.95,OCV_V,linestyle=':',linewidth=lw+0.3,c='k',label='OCV')



#calculate boundary curve
SOC_label='SOC_1'
V_label = 'V_1'

df_boundary_nonzero = df_boundary_data[df_boundary_data[V_label]>0]
SOC_b = df_boundary_nonzero[SOC_label]
V_b = df_boundary_nonzero[V_label]

#need to calculate spline
minV_SOC_spline = splrep(SOC_b,V_b, s=0.0005)
    
#plot boundary
SOC_test = np.linspace(SOC_min_confident,SOC_max_confident,50)
ax.plot(SOC_test, BSpline(*minV_SOC_spline)(SOC_test)-0.002, color='b',
        linewidth=2,label='V(SOC) boundary')
# cutoff_offset = 0.02
# SOC_test = np.linspace(SOC_min_confident,SOC_max_confident,50).reshape(-1,1)
# ax.plot(SOC_test.flatten()+0.03, gpr.predict(SOC_test)-cutoff_offset,color='b',
#         linewidth=2,label='V(SOC) boundary')

plt.legend(edgecolor='w',fontsize=12,loc='lower right')

ax.set_xlim([-0.01,1.01])
# fig.savefig('18July2023_PlatingOnsets_all_higherThreshold_BOL_wBoundaryOCV.svg',bbox_inches = 'tight', dpi= 300)

#%% Explore the plating onsets data for BOL - what can we learn about the datapoints closest to the boundary?
df_simInfo = df_simInfo_original.copy()
ind_i = filterByOnset(df_simInfo, cap_min=0.0, SOC_max=0.95, V_min=3.8, V_max=4.39)

#first find indices of examples that are close to the voltage boundary
df_simInfo_plating = df_simInfo.loc[ind_i,:]#get examples with plating

#get boundary:
SOH_i = 1
SOC_label='SOC_' + str(SOH_i)
V_label = 'V_' + str(SOH_i)

df_boundary_nonzero = df_boundary_data[df_boundary_data[V_label]>0]
SOC_b = df_boundary_nonzero[SOC_label]
V_b = df_boundary_nonzero[V_label]

#calculate spline
minV_SOC_spline = splrep(SOC_b,V_b, s=0.0005)

cmap = matplotlib.cm.get_cmap('coolwarm_r')
norm = matplotlib.colors.Normalize(vmin=0.82,vmax=1.02)
c = cmap(norm(SOH_i))

#calculate the distance from boundary
onsetSOC = df_simInfo_plating['Onset SOC']
df_simInfo_plating['V to boundary'] = df_simInfo_plating['Onset Voltage']-BSpline(*minV_SOC_spline)(onsetSOC)

df_simInfo_plating = df_simInfo_plating.sort_values(by='V to boundary',ascending=True)

#plot distribution of distances from boundary
sns.displot(df_simInfo_plating, x='V to boundary', hue = 'SOH',
            hue_norm=(0.82,1.02), palette = 'coolwarm_r',alpha=1, legend=False,
            height=4, aspect=2, common_norm=False,binwidth=0.005)

#%% Get list of index lists for each region of interest
list_boundaryIndices = []

upperVtoB = [0.022,0.040,0.055]
lowerVtoB = [0,0.0385,0.052]
SOC_cond = df_simInfo_plating['Onset SOC']>0.25

for i in range(len(upperVtoB)):
    V_cond = df_simInfo_plating['V to boundary']<upperVtoB[i]
    V_cond2 = df_simInfo_plating['V to boundary']>lowerVtoB[i]
    ind_closeToBoundary_i = df_simInfo_plating[V_cond & SOC_cond & V_cond2].index
    list_boundaryIndices.append(ind_closeToBoundary_i)

#%% Fig. S4a,f,k. Plot close to boundary points in yellow, others in grey
SOC_max_confident = 0.95
SOC_min_confident = 0.1
plt.rcParams.update({'font.size': 16})
fs=15
lw = 2

for i in range(len(list_boundaryIndices)):
    ind_closeToBoundary = list_boundaryIndices[i]
    
    #create figure
    fig, ax = plt.subplots(nrows=1,ncols=1)
    fig.set_size_inches(6, 3)
    ax.set_ylabel('Voltage (V)', fontsize=fs)
    ax.set_xlabel('State of Charge', fontsize=fs)
    ax.set_ylim([4.03,4.41])
    ax.set_xlim([-0.02,1])
    
    #plot all stars in grey
    ax.scatter(df_simInfo_plating['Onset SOC'],
               df_simInfo_plating['Onset Voltage'],
               marker='*', color=[0.8, 0.8, 0.8], s=50,linewidth=0.5)
    
    #plot close to boundary in yellow
    ax.scatter(df_simInfo_plating['Onset SOC'][ind_closeToBoundary],
               df_simInfo_plating['Onset Voltage'][ind_closeToBoundary],
               marker='*', color=[1, 1, 0], s=100,edgecolors = 'k',linewidth=0.5)

#%% Fig. S4b-e,g-j,l-o. For points close to boundary, plot rate, T, irrevLi, V vs. SOC
graphiteCap = 3.35
plt.rcParams.update({'font.size': 16})
lines = ["-"]
lw = 1
count = 1
ls='-'
for j in range(len(list_boundaryIndices)):
    ind_closeToBoundary = list_boundaryIndices[j]
    fig, ax = plt.subplots(nrows=4,ncols=1,sharex=True,layout='constrained')
    fig.set_size_inches(6, 7)
    
    for i in np.arange(0,len(ind_closeToBoundary),1):  
        df_i = df_outputs[df_outputs['simulation #']==ind_closeToBoundary[i]]
        df_i = df_i[df_i['SOC_full']<df_simInfo['Onset SOC'][ind_closeToBoundary[i]]]
        
    
        ax[3].plot(df_i['SOC_full'],df_i['V'],linewidth=lw,linestyle=ls,label='Ex. '+str(count))
        ax[2].plot(df_i['SOC_full'],df_i['irrevLi_pct']*graphiteCap*1000/100,linewidth=lw,linestyle=ls)
        ax[1].plot(df_i['SOC_full'],df_i['T'],linewidth=lw,linestyle=ls)
        ax[0].plot(df_i['SOC_full'], df_i['rate'],linewidth=lw,linestyle=ls)
        count = count+1
        
    
    ax[0].set_ylabel('C-rate')
    ax[1].set_ylabel('T (' + chr(176) +'C)')
    ax[2].set_ylabel('Li ($\mu$Ah/cm$^2$)')
    ax[2].set_yticks([0.1,0.2,0.3])
    ax[3].set_ylabel('Voltage')
    ax[3].set_xlabel('State of Charge')
    
    ax[0].set_ylim([1.9,8.2])
    ax[0].set_yticks([2,4,6,8])
    ax[1].set_ylim([8,64])
    ax[3].set_ylim([3.56,4.44])
    ax[3].set_yticks([3.6,4,4.4])
    ax[3].set_xlim([-0.01,1.01])


#%% Investigate points closest to boundary for SOH=0.85, study the parameters causing shift
#get SOH aging data
df_outputs = df_outputs_SOH
df_params = df_params_SOH
df_onsets, V_all, SOC_all = utils.parseData_calcOnsets(df_outputs, onset_threshold)
df_simInfo_original = pd.concat([df_params,df_onsets],axis=1)

df_simInfo = df_simInfo_original.copy()
ind_i = filterByOnset(df_simInfo, cap_min=0.0, SOC_max=0.95, V_min=3.8, V_max=4.39)

#first find indices of examples that are close to the voltage boundary
df_simInfo_plating = df_simInfo.loc[ind_i,:]#get examples with plating

#get boundary:
SOH_i = 0.85
SOC_label='SOC_' + str(SOH_i)
V_label = 'V_' + str(SOH_i)

df_boundary_nonzero = df_boundary_data[df_boundary_data[V_label]>0]
SOC_b = df_boundary_nonzero[SOC_label]
V_b = df_boundary_nonzero[V_label]

#calculate spline
minV_SOC_spline = splrep(SOC_b,V_b, s=0.0005)

#calculate the distance from boundary
onsetSOC = df_simInfo_plating['Onset SOC']
df_simInfo_plating['V to boundary'] = df_simInfo_plating['Onset Voltage']-BSpline(*minV_SOC_spline)(onsetSOC)

df_simInfo_plating = df_simInfo_plating.sort_values(by='V to boundary',ascending=True)
SOH_cond = df_simInfo_plating['SOH']==SOH_i

#plot distribution of distances from boundary
sns.displot(df_simInfo_plating[SOH_cond], x='V to boundary', hue = 'SOH',
            hue_norm=(0.82,1.02), palette = 'coolwarm_r',alpha=1, legend=False,
            height=4, aspect=2, common_norm=False,binwidth=0.005)

#%% get index for the region of interest near boundary
upperVtoB = 0.025
lowerVtoB = 0
SOC_cond = df_simInfo_plating['Onset SOC']>0.15


V_cond = df_simInfo_plating['V to boundary']<upperVtoB
V_cond2 = df_simInfo_plating['V to boundary']>lowerVtoB
ind_closeToBoundary = df_simInfo_plating[V_cond & SOC_cond & SOH_cond & V_cond2].index

#%% Fig. S7 Plot close to boundary points in SOH color, others in grey
SOC_max_confident = 0.95
SOC_min_confident = 0.1
plt.rcParams.update({'font.size': 16})
fs=15
lw = 2

#create figure
fig, ax = plt.subplots(nrows=1,ncols=1)
fig.set_size_inches(6, 3)
ax.set_ylabel('Voltage (V)', fontsize=fs)
ax.set_xlabel('State of Charge', fontsize=fs)
ax.set_ylim([4.03,4.41])
ax.set_xlim([-0.02,1])

#plot all stars in grey
ax.scatter(df_simInfo_plating.loc[SOH_cond,'Onset SOC'],
           df_simInfo_plating.loc[SOH_cond,'Onset Voltage'],
           marker='*', color=[0.8, 0.8, 0.8], s=50,linewidth=0.5)

#plot close to boundary in SOH color
cmap = matplotlib.cm.get_cmap('coolwarm_r')
norm = matplotlib.colors.Normalize(vmin=0.82,vmax=1.02)
c = cmap(norm(SOH_i))

ax.scatter(df_simInfo_plating['Onset SOC'][ind_closeToBoundary],
           df_simInfo_plating['Onset Voltage'][ind_closeToBoundary],
           marker='*', color=c, s=100,edgecolors = 'k',linewidth=0.5)

# fig.savefig('19July2023_onsetsNearBoundarySOH_'+ str(SOH_i)+'.svg',bbox_inches = 'tight', dpi= 300)

#%% Table S1. Analyze the indices close to the boundary - what are the average aging parameters?

#The best and worst aging param values correspond to SOH=1, SOH=0.8 respectively

df_agingParams = df_simInfo_plating.loc[ind_closeToBoundary,'an_expansion_um':'ca_LAM']
df_agingStats = df_agingParams.describe().loc['mean':'std',:].transpose()
df_agingStats['best'] = [0, 0, 0, 0, 1, 1, 1, 0, 0, 0]
df_agingStats['worst'] = [8e-6, 8e-6, 0.05, 0.05, 0.5, 0.5, 0.8, 0.08, 0.1, 0.2]
df_agingStats['Norm 0 to 1'] = (df_agingStats['mean']-df_agingStats['worst'])/(df_agingStats['best']-df_agingStats['worst'])
df_agingStats['Norm 0.8 to 1'] = np.round((0.8 + 0.2*df_agingStats['Norm 0 to 1']),decimals=2)
df_agingStats['Std Norm'] = np.round(0.2*df_agingStats['std']/abs(df_agingStats['best']-df_agingStats['worst']),decimals=3)


#%% Fig. 5f Function to visualize onset variability for same boundary intersection
def visualizeVandOnset_sameThreshold(df_simInfo, df_outputs, targetSOC, rangeSOC,
                                     fixed_axes=True,ylim=[3.9,4.4],xlim=[-0.02,1],SOC_error=0):
    SOC_cutoff = df_simInfo['iSOC']+df_simInfo['Cap_cutoff'] #no iSOC added here because cap cutoff acounts for error
    atCutoff = (SOC_cutoff > (targetSOC-rangeSOC)) & (SOC_cutoff < (targetSOC+rangeSOC))
    ind_i_unfiltered = df_simInfo[atCutoff].index
    
    ind_filtered = filterByOnset(df_simInfo,cap_min=0, SOC_max=0.95, V_min=4.0, V_max=4.39,sameRate=True)
    ind_i_filtered = list(i for i in ind_i_unfiltered if i in ind_filtered)
    
    fig, ax = plt.subplots(nrows=1,ncols=1)
    fig.set_size_inches(3, 3)
    plt.rcParams.update({'font.size': 14})
    # ax.plot(OCV_x/max(OCV_x),OCV_V,linestyle=':',linewidth=3,c='k')
    
    #plot voltage curves up to onset 
    for i in ind_i_filtered:
        onset_i = df_simInfo['Onset SOC'][i]+df_simInfo['iSOC_err'][i]
        df_outputs_i = df_outputs[df_outputs['simulation #']==i]
        V_i = df_outputs_i['V']
        SOC_i = df_outputs_i['chargeCap']+df_simInfo['iSOC'][i]+df_simInfo['iSOC_err'][i]
        beforeOnset = SOC_i<onset_i
        
        ax.plot(SOC_i[beforeOnset],V_i[beforeOnset],linewidth=1.2)
            
    ax.scatter(df_simInfo['Onset SOC'][ind_i_filtered]+df_simInfo['iSOC_err'][ind_i_filtered], \
               df_simInfo['Onset Voltage'][ind_i_filtered], marker='*', \
                   color=[1, 1, 0], s=100, edgecolors = 'k',zorder=2 ,label = "Plating Onset") 
    
    SOC_spline = np.linspace(0.02,0.95,50)
    ax.plot(SOC_spline+SOC_error, BSpline(*minV_SOC_spline)(SOC_spline)-y_cutoff_offset, color='b',
            linewidth=2,linestyle='-')
    
    ax.set_ylabel('Voltage (V)')
    ax.set_xlabel('State of Charge')
    if fixed_axes==True:
        ax.set_ylim(ylim)
        ax.set_xlim(xlim)
        #hard code the xticks for paper figure example
        if targetSOC==0.45:
            ax.set_xticks([0.2,0.4,0.6])
            
    else:
        if ind_i_filtered:
            # ax.set_ylim([3.9,max(df_simInfo['Onset Voltage'][ind_i_filtered])+0.05])
            # ax.set_xlim([min(df_simInfo['Onset SOC'][ind_i_filtered])-0.2, \
            #              max(df_simInfo['Onset SOC'][ind_i_filtered])+0.1])
                
            ax.set_ylim([3.9,max(df_simInfo['Onset Voltage'][ind_i_filtered])+0.05])
            ax.set_xlim([min(df_simInfo['Onset SOC'][ind_i_filtered])-0.2, \
                         max(df_simInfo['Onset SOC'][ind_i_filtered])+0.05])
    
    # ax.text
    plt.title('cutoff at ' +str(targetSOC)+' SOC')

    return fig, ax

def addFeatsBeforeCutoff(df_simInfo, df_outputs, nPts, var_name, prefix=""):
    # n_simulations = len(df_simInfo)

    #for each example, add computed features to lists, then add to dataframe at end
    list_var_array_all = []

    # for i in range(n_simulations):
    for i in df_simInfo.index:
        df_outputs_i = df_outputs[df_outputs['simulation #']==i]
        var_i = df_outputs_i[var_name]
        Cap_i = df_outputs_i['chargeCap']        
        Cap_cutoff_i = df_simInfo['Cap_cutoff'][i]
        
        capValsBeforeCutoff = np.linspace(0,Cap_cutoff_i, nPts)
        var_array_i = np.interp(capValsBeforeCutoff, Cap_i, var_i)
        
        #store computed values in list to later add to DataFrame
        list_var_array_all.append(var_array_i)

    
    label = var_name + "_" #default column label prefix
    if prefix:
        label = prefix + "_" #custom provided by user
        
    colNames = utils.createColLabels(label,nPts)
    
    #add columns to df_simInfo with features
    data = np.array(list_var_array_all)
    for i in range(nPts):
        df_simInfo[colNames[i]] = data[:,i]

    return df_simInfo

#%% Visualize plating onsets for trials with same SOC cutoff
df_outputs = df_outputs_BOL
df_params = df_params_BOL

df_onsets, V_all, SOC_all = utils.parseData_calcOnsets(df_outputs, onset_threshold)
df_simInfo_original = pd.concat([df_params,df_onsets],axis=1)
df_simInfo_original['iSOC_err'] = 0 #initialize SOC_err to 0

df_simInfo = df_simInfo_original.copy()
ind_i = filterByOnset(df_simInfo, cap_min=0.0, SOC_max=0.95, V_min=3.8, V_max=4.39)

#need to use only data with plating, otherwise some voltage profiles don't cross V(SOC) boundary
df_simInfo_plating = df_simInfo.loc[ind_i,:]
# df_outputs_plating = df_outputs.loc[ind_i,:]


SOC_label='SOC_1'
V_label = 'V_1'

df_boundary_nonzero = df_boundary_data[df_boundary_data[V_label]>0]
SOC_b = df_boundary_nonzero[SOC_label]
V_b = df_boundary_nonzero[V_label]

#need to calculate spline
minV_SOC_spline = splrep(SOC_b,V_b, s=0.0005)


plt.rcParams.update({'font.size': 15})
fig, ax = plt.subplots(nrows=1,ncols=1)
fig.set_size_inches(6, 4)
SOC_max_confident=0.93
SOC_min_confident=0.12
ax.set_ylabel('Voltage (V)')
ax.set_ylim([3.40,4.44])
# ax.set_ylim([3.9,4.2])
ax.set_xlabel('State of Charge')

ax.scatter(df_simInfo['Onset SOC'][ind_i], df_simInfo['Onset Voltage'][ind_i], marker='*', color=[1, 1, 0], s=100, 
            edgecolors = 'k',linewidth=0.5,zorder=2 ,label = "Plating Onset")      

y_cutoff_offset = 0
x_offset = 0
SOC_test = np.linspace(SOC_min_confident,SOC_max_confident,50)
ax.plot(SOC_test+x_offset, BSpline(*minV_SOC_spline)(SOC_test)-y_cutoff_offset, color='b',
        linewidth=2,label='V(SOC) boundary')

#%% Fig. 5f Plot trials with variable SOC onset given the same cutoff
#add features before and after capacity cutoffs to help with plotting
df_simInfo_plating = addVandCap_boundaryIntersect(df_simInfo_plating, df_outputs,
                                                     minV_SOC_spline, 0,0)
df_simInfo_plating = addFeatsBeforeCutoff(df_simInfo_plating, df_outputs, nPts=5, var_name = 'rate')

# targetSOC_list = np.round(np.linspace(0.3,0.6,7),decimals=3)
targetSOC_list = [0.45]
saveFigs = False
SOC_err = 0.0
for i in targetSOC_list:
    fig,ax = visualizeVandOnset_sameThreshold(df_simInfo_plating, df_outputs, targetSOC=i,
                                              rangeSOC=0.003,fixed_axes=True,ylim=[3.76,4.28],
                                              xlim=[0.1,0.65],SOC_error=SOC_err)
    


