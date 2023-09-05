# -*- coding: utf-8 -*-
"""
This version created on Mon Apr 3, 2023

@author: Zachary M Konz

The goal of this script is to understand how Voltage and Li plating respond to 
changes in variable cell parameters or operating conditions. The analysis herein 
corresponds to Manuscript figures Fig. 2-3, Fig. S2-S3.

"""
#%% Import data and packages
import numpy as np
import matplotlib.pyplot as plt
import warnings
import utils
warnings.filterwarnings("ignore")

OCV_data = np.loadtxt("fullOCVexptCo10extended.CSV", delimiter = ',')# C/10 Gr|NMC532 charging curve
OCV_x = OCV_data[:,0].ravel();
OCV_V = OCV_data[:,1].ravel();

#specify column names in data file, use to create DataFrame from COMSOL txt file 
paramNames = ['iSOC', 'T target',
            'T0','T1','T2','T3','T4','T5','T6','T7',
            't1','t2','t3','t4','t5','t6','t7','t8',
            'ramp1','ramp2','ramp3','ramp4','ramp5','ramp6','ramp7','ramp8',
            'rate1','rate2','rate3','rate4']

agingParamNames = ['an_expansion_um','ca_expansion_um','an_void','ca_void','an_i0_factor',
            'ca_i0_factor','cond_factor','an_LiLoss','an_LAM','ca_LAM','SOH']

outputNames = ['t', 'T', 'V', 'revLi', 'irrevLi',
            'capacity','rate','simulation #']

fileName = "sensitivityAnalysis_simulations.txt"
colNames = paramNames+ agingParamNames + outputNames
df_alldata_all, df_params_all, df_outputs_all = utils.createDfFromTxt(fileName, colNames)

onset_threshold = 0.01 #irreversible Li as percentage of graphite capacity, %
df_onsets_all, V_all, SOC_all = utils.parseData_calcOnsets(df_outputs_all, onset_threshold)

#%% Pick which rate, initial SOC, temperature are the baseline
df = df_params_all.copy()
base1_ind = df[(df['iSOC']==0.1) & (df['rate1']==5) & (df['T0']==308.15)].index[0]
base2_ind = df[(df['iSOC']==0.2) & (df['rate1']==4) & (df['T0']==298.15)].index[0]
base3_ind = df[(df['iSOC']==0.05) & (df['rate1']==7) & (df['T0']==318.15)].index[0]

#NOTE: uncomment one of the three lines below to perform the analysis with different baselines
ind_plot = range(base1_ind,base2_ind) #this will produce Fig. 2-3 in manuscript
# ind_plot = range(base2_ind,base3_ind) #this will produce Fig. S2 in Supporting Info
# ind_plot = range(base3_ind,len(df_params_all)) #this will produce Fig. S3 in Supporting Info

#the following conditional sequence changes the plot text label according to baseline
if ind_plot[0]==base1_ind:
    plt_txt_label = '5C, 10% iSOC, 35'+chr(176)+'C'    
elif ind_plot[0]==base2_ind:
    plt_txt_label = '4C, 20% iSOC, 25'+chr(176)+'C'
else:
    plt_txt_label = '7C, 5% iSOC, 45'+chr(176)+'C'
    
df_params = df_params_all.loc[ind_plot,:]
df_onsets = df_onsets_all.loc[ind_plot,:]
df_alldata = df_alldata_all[df_alldata_all['simulation #'].isin(ind_plot)]
df_outputs = df_outputs_all[df_outputs_all['simulation #'].isin(ind_plot)]

#%% Fig. 2a-l. (also Fig. S2,S3) Data visualization: for each parameter, plot the voltage profile for each condition
list_plotTitles = ['T','C-rate','iSOC','Gr expansion','NMC expansion',
                    'Gr voids','NMC voids','Gr exchange current','NMC exchange current',
                    'conductivity','Gr Li loss','Gr loss of active matl','NMC loss of active matl']

list_paramNames = ['T0','rate1','iSOC','an_expansion_um','ca_expansion_um','an_void',
                   'ca_void','an_i0_factor','ca_i0_factor','cond_factor','an_LiLoss','an_LAM','ca_LAM']
            
list_units = [chr(176)+'C','C','% SOC','um','um','%','%','*i0','*i0','*k','%','%','%']

list_multipliers = [1,1,100,1000000,1000000,100,100,1,1,1,100,100,100]

baseline_ind = df_outputs['simulation #']==df_outputs['simulation #'].iloc[0] #the first simulation is the baseline
V_baseline = df_outputs[baseline_ind]['V']
SOC_baseline = df_outputs[baseline_ind]['SOC_full']

plotAsSubplots = True
saveFigure = False
saveIndividualFigures = False
plotChargeCapacity = True

if plotAsSubplots:
    n_wide = 3;
    # n_high = 5;
    n_high = 4;
    fs = 12;
    fig, axs = plt.subplots(nrows=n_high,ncols=n_wide,sharex='col',sharey='row',layout='constrained')
    # fig.set_size_inches(5*n_wide,3*n_high)
    fig.set_size_inches(4.3*n_wide,2.5*n_high)

list_num = range(0,len(list_paramNames))
leave_out = 4 #this corresponds to NMC electrode expansion
ind_toPlot = list(i for i in list_num if i!=leave_out)
count = 0
# for i in range(0,len(list_paramNames)):
for i in ind_toPlot:
    if plotAsSubplots:
        ax_i = axs.flat[count]
        count = count + 1
    else:
        fig, ax_i = plt.subplots(nrows=1,ncols=1,sharex='col',sharey='row',layout='constrained')
        fig.set_size_inches(5,3)
        fs = 12

    
    #plot 1st simulation with standard values
    SOC_base_plot = SOC_baseline
    
    if plotChargeCapacity:
        SOC_base_plot = SOC_base_plot-SOC_baseline.iloc[0]
        
    ax_i.plot(SOC_base_plot, V_baseline, ':k',label='baseline')
    ax_i.set_prop_cycle(None)
    
    #find simulation numbers corresponding to non-standard values
    param_i = df_params[list_paramNames[i]]
    std_value_i = param_i.iloc[0]
    nonStd = param_i!=std_value_i
    simulationNums_i = df_params[nonStd].index.tolist()
    
    for j in simulationNums_i:
        data_j = df_outputs[df_outputs['simulation #']==j]
        param_value = param_i[j]
        display_val = list_multipliers[i]*param_value
        initial_SOC = data_j['SOC_full'].iloc[0]
        SOC_plot = data_j['SOC_full']

        if plotChargeCapacity:
            SOC_plot = SOC_plot-initial_SOC
        
        #reduce number of decimals in label
        if display_val >= 1:
            display_val = int(display_val)
        
        if count==1:#if temperature, then convert to degrees Celsius
            label_i = str(display_val-273)+list_units[i]
        else:
            label_i = str(display_val)+list_units[i]
        
        ax_i.plot(SOC_plot,data_j['V'],label=label_i)
    
    ax_i.set_ylim([3.5,4.39])
    ax_i.set_xlim([0,0.79])
    ax_i.text(0.02, 4.3, list_plotTitles[i],fontsize = fs)
    
    x_axis_label = 'SOC'
    if plotChargeCapacity:
        x_axis_label = 'Charge capacity, SOC-SOC$_0$'
    
    #plot formatting based on whether subplots or individual figures
    if plotAsSubplots:
        if (count+2)%3 == 0: #only label y axis for left column of plots
            ax_i.set_ylabel('Gr|NMC Voltage (V)',fontsize=fs)
        if count > 9: #only label x axis for bottom column of plots
            ax_i.set_xlabel(x_axis_label,fontsize=fs)
    else:
        ax_i.set_ylabel('Gr|NMC Voltage (V)',fontsize=fs)
        ax_i.set_xlabel(x_axis_label,fontsize=fs)
    
    ax_i.legend(loc='lower right', fontsize= (fs-2),edgecolor='w')
    
    #plot onsets for each condition after adding curves
    #plot baseline condition
    initial_SOC_base = SOC_baseline.iloc[0]
    onset_SOC_base = df_onsets['Onset SOC'].iloc[0]
    
    if plotChargeCapacity:
        onset_SOC_base = onset_SOC_base - initial_SOC_base
    
    ax_i.scatter(onset_SOC_base, df_onsets['Onset Voltage'].iloc[0], 
                 marker='*',  c='k',s=200, zorder=2, label=None) 

    ax_i.set_prop_cycle(None)
    for j in simulationNums_i:
        data_j = df_onsets.loc[j,:]
        param_value = param_i[j]
        initial_SOC = df_outputs[df_outputs['simulation #']==j]['SOC_full'].iloc[0]
        onset_SOC = data_j['Onset SOC']
        
        if plotChargeCapacity:
            onset_SOC = onset_SOC-initial_SOC

        ax_i.scatter(onset_SOC, data_j['Onset Voltage'], marker='*', s=200, zorder=2, label=None)
    
    x_axis = '_vsSOC'
    
    if plotChargeCapacity:
        x_axis = '_vsChargeCap'
        
    if saveIndividualFigures:
        fileName = "agingParamsTest_updatedEqns_5Cbaseline_19June2023" + x_axis + list_paramNames[i] +".svg"
        plt.savefig(fileName,bbox_inches = 'tight', dpi= 300) 

if saveFigure:
    fileName = "sensitivityAnalysis" + str(plt_txt_label) + "_editDate.svg"
    # plt.savefig(fileName,bbox_inches = 'tight', dpi= 300)    

#%% Functions for calculating voltage differences
from scipy import interpolate

def calculateDiff_VQ(SOC1, V1, SOC2, V2):
    """
    Input: 2 voltage-SOC curves for charging
      SOC1 (pd.Series): full-cell SOC for charging condition 1
      V1 (pd.Series): full-cell voltage for charging condition 1
      SOC2 (pd.Series): full-cell SOC for charging condition 2
      V2 (pd.Series): full-cell voltage for charging condition 2

    Output:
      np.mean(diff_xtoy): mean voltage difference condition1-condition2 between x-y% capacity
      np.mean(diff_xxxtoy): mean capacity difference condition1-condition2 between xxx-y voltage
    """    
    cap1 = SOC1-min(SOC1)
    cap2 = SOC2-min(SOC2)
    
    nPts = 50
    Cap_range = [0,0.02]
    Q_Oto2 = np.linspace(Cap_range[0],Cap_range[1],nPts)
    V1_interp_0to2 = interpolate.splev(Q_Oto2, interpolate.splrep(cap1,V1), der=0)
    V2_interp_0to2 = interpolate.splev(Q_Oto2, interpolate.splrep(cap2,V2), der=0)
    diff_0to2 = V1_interp_0to2-V2_interp_0to2
    
    Cap_range = [0.02,0.05]
    Q_2to5 = np.linspace(Cap_range[0],Cap_range[1],nPts)
    V1_interp_2to5 = interpolate.splev(Q_2to5, interpolate.splrep(cap1,V1), der=0)
    V2_interp_2to5 = interpolate.splev(Q_2to5, interpolate.splrep(cap2,V2), der=0)
    diff_2to5 = V1_interp_2to5-V2_interp_2to5
    
    V_range = [3.95,4]
    V_395to4 = np.linspace(V_range[0],V_range[1],nPts)
    cap1_interp_395to4 = interpolate.splev(V_395to4, interpolate.splrep(V1,cap1), der=0)
    cap2_interp_395to4 = interpolate.splev(V_395to4, interpolate.splrep(V2,cap2), der=0)
    diff_395to4 = cap1_interp_395to4-cap2_interp_395to4
    
    return np.mean(diff_0to2),np.mean(diff_2to5),np.mean(diff_395to4)

def getVandOnsetChanges(simNums,df_outputs,df_onsets):
    """
    For 2 simulations that vary only by one parameter change,
    calculate the voltage curve shift and plating onset change 
    
    Input: 2 voltage-SOC curves for charging
      simNums (tuple): two unique simulation numbers [n1, n2]
      df_outputs (DataFrame): 
      df_onsets (DataFrame):

    Output:
      d1: mean voltage difference n1-n2 between 0-2% capacity
      d2: mean voltage difference n1-n2 between 2-5% capacity
      d3: mean capacity difference n1-n2 between 3.95-4.00 V
      diff_SOC: onset SOC difference n1-n2
      diff_V: onset voltage difference n1-n2
      
    """ 
    n1 = simNums[0]
    n2 = simNums[1]
    
    #get voltage-SOC curve data for both parameter values
    SOC1 = df_outputs[df_outputs['simulation #']==n1]['SOC_full']
    V1 = df_outputs[df_outputs['simulation #']==n1]['V']

    SOC2 = df_outputs[df_outputs['simulation #']==n2]['SOC_full']
    V2 = df_outputs[df_outputs['simulation #']==n2]['V']
    
    #feed into curves data
    d1, d2, d3 = calculateDiff_VQ(SOC1, V1, SOC2, V2)#note, differences are calculated n1-n2, order matters
    diff_SOC = df_onsets.loc[n1,'Onset SOC']-df_onsets.loc[n2,'Onset SOC']
    diff_V = df_onsets.loc[n1,'Onset Voltage']-df_onsets.loc[n2,'Onset Voltage']
    
    return d1,d2,d3,diff_SOC,diff_V
    
#%% Get voltage shift for all parameter change combinations within delta (new code 6/5/23)
param_delta_max = [10,1,0.05,4e-6,4e-6,0.03,0.03,0.4,0.4,0.2,0.05,0.06,0.15]

import itertools

#nested lists, length=n_params, each element has variable length=n_combinations
list_diff_0to2 = []
list_diff_2to5 = []
list_diff_395to4 = []
list_onset_diff_V = []
list_onset_diff_SOC = []
list_simNums = []#each element is a list of tuples corresponding to simulation numbers

for i in range(len(param_delta_max)):
    #find simulation numbers corresponding to changing only one parameter
    index_baseline = ind_plot[0]#scalar
    param_i = df_params[list_paramNames[i]]
    std_value_i = param_i[index_baseline]
    index_changes = df_params[param_i!=std_value_i].index.tolist()
    
    index_list_i = [index_baseline]+index_changes #list of index of simulations of interest
    
    #find all unique pairs of simulations that meet criteria diff(param1,param2) < param_delta_max
    list_pairs = list(itertools.combinations(index_list_i,2))
    list_pairs_filtered_i = []
    
    for k in range(len(list_pairs)):
        value1 = param_i[list_pairs[k][0]]
        value2 = param_i[list_pairs[k][1]]
        
        if abs(value1-value2)<=param_delta_max[i]:
            list_pairs_filtered_i.append(list_pairs[k])

    #for each pair, calculate onset and voltage differences, add to lists
    list_d0to2_i = []
    list_d2to5_i = []
    list_d395to4_i = []
    list_SOCdiff_i = []
    list_Vdiff_i = []
    
    for j in range(len(list_pairs_filtered_i)):
        simNums = list_pairs_filtered_i[j]#tuple [n1,n2]
        d1, d2, d3, onsetSOC_shift, onsetV_shift = getVandOnsetChanges(simNums,df_outputs,df_onsets)
        
        list_d0to2_i.append(d1)
        list_d2to5_i.append(d2)
        list_d395to4_i.append(d3)
        list_SOCdiff_i.append(onsetSOC_shift)
        list_Vdiff_i.append(onsetV_shift)
    
    list_diff_0to2.append(list_d0to2_i)
    list_diff_2to5.append(list_d2to5_i)
    list_diff_395to4.append(list_d395to4_i)
    list_onset_diff_V.append(list_Vdiff_i)
    list_onset_diff_SOC.append(list_SOCdiff_i)
    list_simNums.append(list_pairs_filtered_i)

#%% Fig. 3c-f, (also Fig. S2,S3) Plot voltage shifts all parameter combinations
list_paramLabels = ['Temp.','Rate','iSOC','Gr expand','NMC expand',
                    'Gr void','NMC void','Gr i0','NMC i0',
                    'Conductivity','Gr Li loss','Gr LAM','NMC LAM']

plotAsSubplots = True
saveIndividualFigures = False
plotZoom = False

pt_size = 12
if plotAsSubplots:
    # n_wide = 3;
    n_wide = 2;
    n_high = 2;
    fs = 12.5;
    fig, axs = plt.subplots(nrows=n_high,ncols=n_wide,sharex='col',sharey='row',layout='constrained')
    fig.set_size_inches(3*n_wide,3*n_high)

list_metrics = [list_diff_2to5,list_diff_395to4]
metric_labels = ['Voltage shift, $\Delta$V$_{2-5\%}$ (V)','Voltage shift, $\Delta$Q$_{4V}$ (SOC)']
x_limits_metrics = [[0,0.068],[0,0.13]]
x_limits_metrics_zoom = [[0,0.02],[0,0.04]]
x_onset_ticks_zoom = [[0,0.02,0.04],[0,0.02,0.04]]

list_onsets = [list_onset_diff_V,list_onset_diff_SOC]
onset_labels = ['$\Delta$Onset, Voltage (V)','$\Delta$Onset, Capacity (SOC)']
y_limits_onsets = [[-0.11,0.11],[-0.2,0.2]]
y_limits_onsets_zoom = [[-0.03,0.03],[-0.03,0.08]]
y_onset_ticks = [[-0.1,-0.05,0,0.05,0.1],[-0.20,-0.1,0,0.1,0.2]]
y_onset_ticks_zoom = [[-0.02,0,0.02],[0,0.05]]
y_var = []

for i in range(len(list_metrics)):
# for i in [1,2]:#only plot voltage response 2-5, 3.95-4.0V 
    for j in range(len(list_onsets)):
        if plotAsSubplots:
            # ax_i = axs[i,j]
            ax_i = axs[j,i]
        else:
            fig, ax_i = plt.subplots(nrows=1,ncols=1,sharex='col',sharey='row',layout='constrained')
            if plotZoom:
                fig.set_size_inches(2.4,2) 
                pt_size = 18
            else:
                fig.set_size_inches(4,4)
            fs = 12

        #Plot y=0 line on all plots to give visual reference
        ax_i.plot([-1,1],[0,0],'k:',linewidth=1, label=None, zorder=0)
        
        #For each plot, iterate through all parameters, combinations and plot

        marker_list = ['o', '^', '*', 's', 'h', '>','d','p', 'X', 'P', 'D', 'v','2']
        
        cmap = plt.cm.get_cmap('tab10')
        color_list = cmap(np.linspace(0,1,len(marker_list)))
        np.random.seed(seed=1)
        np.random.shuffle(color_list)
        
        for k in range(len(list_paramNames)):
            for l in range(len(list_metrics[0][k])):#for all combinations of parameters
                diff_metric = list_metrics[i][k][l]
                onset_metric = list_onsets[j][k][l]
                
                #check signs of voltage and onset change; if opposite signs, blue, else, red
                
                if diff_metric*onset_metric<0:
                    plotColor = 'b'
                    onset_metric = -1*abs(onset_metric)
                else:
                    plotColor = 'r'
                    onset_metric = abs(onset_metric)
                    
                if abs(onset_metric)<0.001:
                    plotColor = 'gray'
                
                if l==0:
                    lab = list_paramLabels[k]
                else:
                    lab = None
                    
                ax_i.scatter(abs(diff_metric), onset_metric , marker=marker_list[k],
                             s=pt_size, c=[color_list[k]], label=lab)
            
        #Plot formatting depends on whether plotting as subplots or individual
        if plotAsSubplots:
            #only label far left y axis and bottom x axis
            if j==1:
                ax_i.set_xlabel(metric_labels[i], fontsize=fs)
            
            if i==0:
                ax_i.set_ylabel(onset_labels[j], fontsize=fs)
                ax_i.set_yticks(y_onset_ticks[j])
            
            #set x,y limits the same for all plots
            ax_i.set_xlim(x_limits_metrics[i])
            ax_i.set_ylim(y_limits_onsets[j])
        
            if (i==0 & j==0):#label baseline conditions in upper left plot
                ax_i.text(.01, .99, str(plt_txt_label), ha='left', va='top', transform=ax_i.transAxes, fontsize=fs-2)
            
        else: #if plotting indiviual plots
            if plotZoom:
                fs = 9
                
            #labels axes of all plots
            ax_i.set_ylabel(onset_labels[j], fontsize=fs)
            ax_i.set_xlabel(metric_labels[i], fontsize=fs)
            
            if plotZoom:
                ax_i.set_xlim(x_limits_metrics_zoom[i])
                ax_i.set_ylim(y_limits_onsets_zoom[j])
                ax_i.set_yticks(y_onset_ticks_zoom[j])
                ax_i.set_xticks(x_onset_ticks_zoom[i])
                
            else:
                ax_i.set_xlim(x_limits_metrics[i])
                ax_i.set_ylim(y_limits_onsets[j])
                ax_i.set_yticks(y_onset_ticks[j])
            
        if saveIndividualFigures:
            fileName = '9July2023' + str(plt_txt_label) + str(i)+str(j)+".svg"
            plt.savefig(fileName,bbox_inches = 'tight', dpi= 300)
  
        
# fileName = 'SensitivityAnalysis.svg'
# plt.savefig(fileName,bbox_inches = 'tight', dpi= 300)

#%% Fig. 3a. Plot 2 voltage curves for illustration figure

fig, ax_i = plt.subplots(nrows=1,ncols=1,sharex='col',sharey='row',layout='constrained')
fig.set_size_inches(3,3)
fs = 12

    
#find simulation numbers corresponding to desired conditions
#want 20degC, 40degC for baseline 5C, 10% iSOC data

sim_ind1 = df[(df['iSOC']==0.1) & (df['rate1']==5) & (df['T0']==303.15)].index[0]
sim_ind2 = df[(df['iSOC']==0.1) & (df['rate1']==5) & (df['T0']==313.15)].index[0]
ind_plot = [sim_ind1,sim_ind2]

df_params = df_params_all.loc[ind_plot,:]
df_onsets = df_onsets_all.loc[ind_plot,:]
df_alldata = df_alldata_all[df_alldata_all['simulation #'].isin(ind_plot)]
df_outputs = df_outputs_all[df_outputs_all['simulation #'].isin(ind_plot)]

simulationNums_i = ind_plot

for j in simulationNums_i:
    data_j = df_outputs[df_outputs['simulation #']==j]
    param_value = param_i[j]
    display_val = list_multipliers[i]*param_value
    initial_SOC = data_j['SOC_full'].iloc[0]
    SOC_plot = data_j['SOC_full']
    SOC_plot = SOC_plot-initial_SOC
    
    ax_i.plot(SOC_plot,data_j['V'],label=None, c='k',linewidth=1)

ax_i.set_ylim([3.6,4.3])
ax_i.set_xlim([0,0.6])
ax_i.set_yticks([3.6,3.8,4,4.2])
ax_i.set_ylabel('Gr|NMC Voltage (V)',fontsize=fs)
    
x_axis_label = 'Charge capacity, SOC-SOC$_0$'

ax_i.set_xlabel(x_axis_label,fontsize=fs)

#plot onsets for each condition after adding curves

ax_i.set_prop_cycle(None)
for j in simulationNums_i:
    data_j = df_onsets.loc[j,:]
    param_value = param_i[j]
    initial_SOC = df_outputs[df_outputs['simulation #']==j]['SOC_full'].iloc[0]
    onset_SOC = data_j['Onset SOC']
    onset_SOC = onset_SOC-initial_SOC

    ax_i.scatter(onset_SOC, data_j['Onset Voltage'], c='k',marker='*', s=100, zorder=2, label=None)

fileName = 'Vshift_illustration.svg'
plt.savefig(fileName,bbox_inches = 'tight', dpi= 300) 
