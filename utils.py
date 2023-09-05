#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 10 2023

@author: Zachary M Konz
"""
import pandas as pd
import numpy as np

#%% constants:
D = 1.4 
A_c = np.pi*D**2/4
cathodeCap = 2.80 # mAh/cm2
graphiteCap = 3.35 # mAh/cm2
Li_umol_to_mAh = 0.0374 #mAh per umol Li

#%% functions
def createDfFromTxt(fileName, colNames, iSOC=0):
    all_data = np.loadtxt(fileName,skiprows=5)
    tCol = colNames.index('t')
    
    simulation_params,idx = np.unique(all_data[:,:tCol], axis=0, return_index=True)
    simulation_params = simulation_params[np.argsort(idx)]
    n_simulations = np.shape(simulation_params)[0]
    
    #add a column that indexes each simulation from 0 to (n_simulations-1)
    all_data = np.concatenate([all_data,np.zeros((len(all_data),1))],axis=1)
    idxCol = all_data.shape[1]-1
    
    for i in range(n_simulations):
        params_i = simulation_params[i,:]
        idx = np.where((all_data[:,:tCol] == params_i).all(axis=1))
        all_data[idx,idxCol] = i
    
    df_alldata = pd.DataFrame(data=all_data, columns = colNames)
    df_params = pd.DataFrame(data=simulation_params,columns=colNames[:tCol])
    
    if 'iSOC' in colNames: 
        df_outputs = df_alldata.copy()[['simulation #','iSOC','t','rate','capacity','V','T']]
        #Add column for SOC, defined by limiting cathode (full-cell SOC)
        df_outputs['SOC_full'] = df_alldata['iSOC'] + df_alldata['capacity']/A_c/cathodeCap
    else:
        df_outputs = df_alldata.copy()[['simulation #','t','rate','capacity','V','T']]
        #Add column for SOC, defined by limiting cathode (full-cell SOC)
        df_outputs['SOC_full'] = iSOC + df_alldata['capacity']/A_c/cathodeCap
    
    df_outputs['irrevLi_pct'] = df_alldata['irrevLi']*Li_umol_to_mAh*100/graphiteCap
    df_outputs['chargeCap'] = df_alldata['capacity']/A_c/cathodeCap
    
    return df_alldata, df_params, df_outputs

def parseData_calcOnsets(data, threshold):    
    """
    Input:
      data (DataFrame):  with columns 'simulation #', 'SOC_full', 'V', 'irrevLi_pct' 
      threshold (float): irrev. Li threshold used to define plating onset, % of graphite capacity

    Output:
      onsetV (1D array): onset Voltages for all simulations
      onsetSOC (1D array): onset capacity, or full-cell SOC at onset
      V_all (list of 1D arrays): V(t) for each simulation
      SOC_all (list of 1D arrays): full-cell SOC(t) for each simulation
    """
    
    #initialize function output variables
    n_simulations = int(data['simulation #'].max()+1)
    # n_params = np.shape(params)[1]
    V_all = []
    SOC_all = []
    onsetV = np.zeros((n_simulations,1))
    onsetSOC = np.zeros((n_simulations,1))
    onsetCap = np.zeros((n_simulations,1))
    onsetRate = np.zeros((n_simulations,1))
    sim_num = np.zeros(n_simulations)
    
    for i in range(0,n_simulations):
        # indices = np.where((data[:,0:n_params] == params[i,:]).all(axis=1))
        indices = data['simulation #']==i
        data_i = data[indices]
        V_i = data_i['V'].to_numpy()
        Q_i = data_i['SOC_full'].to_numpy()
        irrevLi_i = data_i['irrevLi_pct'].to_numpy()
        rate_i = data_i['rate'].to_numpy()
        V_all.append(V_i)
        SOC_all.append(Q_i)
        
        #calculate plating onset Q, V where irrevLi_i equals threshold by interpolation of simulation data
        if max(irrevLi_i)>threshold:
            onsetV[i] = np.interp(threshold, irrevLi_i.ravel(), V_i.ravel()) # need to use .ravel() to make 1D sequence for fn input
            onsetSOC[i] = np.interp(threshold, irrevLi_i.ravel(), Q_i.ravel())
            onsetCap[i] = onsetSOC[i]-Q_i[0]
            onsetRate[i] = np.interp(threshold, irrevLi_i.ravel(), rate_i.ravel())
        else:
            onsetV[i] = 4.4
            onsetSOC[i] = 1.1
            onsetCap[i] = 1
            onsetRate[i] = 0
            
        sim_num[i] = i
        
    #Create dataframe with onsets information only, useful for EDA
    onsetsData = np.concatenate((onsetV,onsetSOC,onsetCap,onsetRate),axis=1)

    df_onsets = pd.DataFrame(data = onsetsData, index= sim_num.astype(int).tolist(),
                             columns = ['Onset Voltage','Onset SOC','Onset Cap','Onset rate'])
    
    return df_onsets, V_all, SOC_all

def createColLabels(prefix, n):
    colNames = []
    for k in range(n):
        colNames.append(prefix+str(k))
        
    return colNames
