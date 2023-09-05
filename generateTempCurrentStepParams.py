# -*- coding: utf-8 -*-
"""
Generating joint Temperature and Current demand profiles for Charging
Author: Zachary M. Konz

The goal of this code is to generate diverse, realistic, Li-plating-inducing
fast charge protocols (I and T demand) as inputs for a P2D COMSOL simulation. 
This code was used to generate the charge protocols analyzed in Manuscript
Figs. 4-6 and S4-S7. A written description of this code is provided in the 
Supplementary Information.

We have defined the I(t) and T(t) as piecewise functions in the COMSOL software
with 4 and (4*2=) 8 segments respectively. Slight edits could generalize this code
for any number of current or temperature steps.

For some simulations, aging parameters are inlcuded. The aging parameter values
are correlated to the overall battery health, or State-of-Health (SOH).
This value ranges 0.8 to 1.0, and is referred to as 'health index' in the Manuscript.
"""
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import random
from matplotlib.lines import Line2D
import pandas as pd

#%% Specify simulation parameters and script options


rates_max = [8, 7, 6, 5] # maximum allowable current for steps 1, 2, 3, 4 respectively
rates_min = [3, 3, 2, 2] # minimum current

includeAgingParams = True  # if True, include aging parameters in output
agingCorrelatedSOH = False # if True, aging parameters will be correlated relative to SOH/health-index
noAging = True #if True, include default aging parameter values, SOH/health-index=1

agingParamNames = ['an_len_expansion', 'ca_len_expansion', 'an_void_eps', 'ca_void_eps',
                   'an_i0_factor', 'ca_i0_factor', 'el_kappa_factor', 'an_LiLoss',
                   'an_LAM', 'ca_LAM']

agingParamBestVal = [0, 0, 0, 0, 1, 1, 1, 0, 0, 0]
agingParamWorstVal = [8e-6, 8e-6, 0.05, 0.05, 0.5, 0.5, 0.8, 0.08, 0.1, 0.2]

# for correlated aging parameters
SOH_values = [0.98, 0.95, 0.90, 0.85]
SOH_min = 0.8
SOH_max = 1.0

# the rate of the 4th current step should not exceed the 3rd by more than maxRate4_increase
# this constraint makes protocol more realistic; no big current increase near end of charge
maxRate4_increase = 0.5

n_steps = len(rates_min)  # the number of current steps
n_split = 2  # number of temperature steps per each current step
n_temp_steps = n_split*n_steps  # number of temperature steps

T_ramp_ref_max = 20 # max temperature ramp rate in degrees C per minute, at rate=8C
T_ramp_ref_min = 5  # min ramp rate
ref_rate = 8  # reference C-rate for these ramps, used to scale T ramps according to present rate

T_initial_min = 10  # min initial charging temperature
T_initial_max = 45  # max "
T_target_max = 60 # maximum target temperature - T will approach this value, and if reached, fluctuate nearby
T_target_min = 30  # minimum "

# the minimum temperature difference between initial and target temperatures
# alternatively, this means that T_target will always be >= T_initial + T_delta_min
T_delta_min = 5

# T drift rate fraction of the max T ramp rate
# After T_target is reached, the T should fluctuate proportional to applied current
T_drift_factor_max = 0.2

# if it is desired to run simulations only at discrete SOC values, they can be selected from below list
discreteSOC = False
if discreteSOC:
    discrete_SOC_vals = [0.02, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4]
else:  # specify range of initial SOC, which will be randomly selected, rounded to 0.01
    iSOC_min = 0.02
    iSOC_max = 0.5

stopSOC = 0.95 #charge target stop SOC

G = 1/60  # SOC per minute for 1C rate
upper_thresh = 2  # deg C
lower_thresh = 1  # deg C

#%% Functions for generating T, Current profiles

def getInitialRamp(T_initial, rate_initial):
    """
    Goal: select initial T ramp (degC/min) for protocol given the initial T and rate.
    For lower initial T, ramps should generally be larger (rapid heating) to reach high T
    typically desired for fast charge performance. Possible ramp values should be higher
    for higher rates, roughly scaling with I^2*R

    Input:
      T_initial (float): random initial temperature for given simulation
      rate_initial (float): random initial rate for given simulation

    Output:
      ramp_initial (float): initial T ramp rate during charge, degrees C/min

    Constants defined in 1st cell:
      T_initial_min (float): minimum bound of T parameter range
      T_initial_max (float): maximum bound of T parameter range
      T_ramp_ref_min (float): minimum bound of T ramp at reference rate, degrees C/min
      T_ramp_ref_max (float): maximum bound of T ramp at reference rate, degrees C/min
    """

    # Normalize initial T between 0 and 1 according to min and max values
    T_init_norm = (T_initial-T_initial_min)/(T_initial_max-T_initial_min)

    # Select initial ramp, normalized between 0 and 1
    mean = 1-T_init_norm  # randomly selected ramp rate should be higher if T_initial is lower
    sd = 0.4
    ramp_init_norm = np.random.normal(mean, sd, 1)[0]

    if (ramp_init_norm > 1):
        ramp_init_norm = 1
    elif (ramp_init_norm < 0):
        ramp_init_norm = 0

    # scale the ramp rate to the given initial c-rate, scale by ratio of I^2/Iref^2
    ramp_max = (rate_initial**2)/(ref_rate**2)*T_ramp_ref_max
    ramp_min = (rate_initial**2)/(ref_rate**2)*T_ramp_ref_min
    ramp_initial = ramp_min + ramp_init_norm*(ramp_max-ramp_min)

    return np.round(ramp_initial, decimals=1)

def getRandomDrift(rate):
    """
    Goal: after the target T is reached, a random T ramp (degC/min) should be selected
    to mimic random T fluctuation due to imperfect battery thermal management system.
    The possible ramp range should larger for more rapid charging (higher C-rates)

    Input:
      rate (float): current rate applied at this specific step

    Output:
      drift (float): temperature drift rate during charge, degrees C/min

    Constants defined in 1st cell:
      T_drift_factor_max (float): maximum fraction of ramp rate designated for temperature drift
      T_ramp_ref_max (float): maximum bound of T ramp at reference rate, degrees C/min

    """
    mean = 0  # ramp can be positive or negative - negative is possible with battery cooling
    sd = T_drift_factor_max/3  # to get 99% of drift rates within maximum drift factor

    drift_norm = np.random.normal(mean, sd, 1)[0]

    ramp_max = (rate**2)/(ref_rate**2) * T_ramp_ref_max  # maximum expected ramp rate

    drift = drift_norm*ramp_max  # drift is a fraction of the maximum expected ramp

    return drift

def getSOCsteps(n_steps, SOC_initial, stopSOC):
    # get random SOC length for each step, adding variability to multi-step CC
    # randomly generate 8 numbers in same range, then divide by sum
    randNum = np.random.uniform(0.02, 1, n_steps)
    randNum_norm = randNum/sum(randNum)  # now they all add to 1

    return np.round(randNum_norm*(stopSOC-SOC_initial), decimals=3)


#%% Generate profiles
n_conditions = 1000  # number of parameter sets to generate

# initialize lists to store variable lists, will have length of n_conditions
list_rates = []
list_temps = []
list_ramps = []
list_ts = []
list_SOCs = []
list_T_target = []

for j in range(n_conditions):

    # select 4 current step values for jth simulation
    rates = []

    for i in range(len(rates_max)):  # for each current step

        # the maximum possible rate for each step is previously defined in rates_max
        # for the last (4th) step, the rate should not be much larger than previous value
        if i < (len(rates_max)-1):
            rates_max_i = rates_max[i]
        else:
            rates_max_i = min([(rates[i-1]+maxRate4_increase), rates_max[i]])

        # use uniform distribution to select a random rate for each step
        rangeRates = rates_max_i - rates_min[i]
        rand_rate = rates_min[i] + random.uniform(0, 1)*rangeRates
        rates.append(np.round(rand_rate, decimals=1))

    # select initial SOC value
    if discreteSOC:
        SOC_initial = random.choice(discrete_SOC_vals)
    else:
        SOC_initial = np.round(np.random.uniform(
            iSOC_min, iSOC_max), decimals=2)

    # select an initial charging temperature:
    T_initial = np.round(T_initial_min + random.uniform(0, 1)
                         * (T_initial_max-T_initial_min))

    # based on the initial charging temperature and current rate, get initial T ramp rate
    ramp_initial = getInitialRamp(T_initial, rates[0])

    # Next step: pick the target temperature. The target must be higher than the
    # initial temperature by at least T_delta_min

    if T_initial < (T_target_min - T_delta_min):
        # Pick random target temperature between min and max specified target T
        T_target = T_target_min + \
            random.uniform(0, 1)*(T_target_max-T_target_min)
    else:
        # Pick random target temperature between (T_initial+ T_delta_min) and T_target_max
        T_target = (T_initial+T_delta_min) + random.uniform(0, 1) * \
            (T_target_max-(T_initial+T_delta_min))

    # randomly determine SOC passed (proportional to charge time) for each half current step
    SOC_deltas = getSOCsteps(n_steps*n_split, SOC_initial, stopSOC)

    # Given all parameters, calculate temperature profiles for each half current step

    # initialize the following variables at the start of charge, will append values to
    # these lists for subsequent steps

    t_i = [0]  # initial time for each step, in minutes
    T_i = [T_initial]
    ramp_i = [ramp_initial]
    # this duplicates each rate in the list, in order
    rate_i = [i for i in rates for r in range(n_split)]
    rate_i.append(rate_i[-1])
    SOC_i = list((SOC_initial+sum(SOC_deltas[:i]))
                 for i in range(n_steps*n_split))
    SOC_i = SOC_i + [stopSOC]  # add last SOC to end

    rate_to_t = SOC_deltas/G  # minutes per 1C rate for each step
    t_delta_i = (rate_to_t/np.array(rate_i[:(n_steps*n_split)])).tolist()
    # duplicate last element to keep i+1 indexing happy in loop
    t_delta_i.append(t_delta_i[-1])

    T_target_reached = False
    for i in range(n_temp_steps):

        # calculate and append time at the end of ith step
        t_i.append(t_i[i]+t_delta_i[i])
        T_delta = ramp_i[i]*t_delta_i[i]  # T change for ith step
        T_final = T_i[i] + T_delta  # T at end of ith step

        # the final temperature of step i is initial temperature of (i+1)
        T_i.append(T_final)

        if T_final > T_target:
            T_target_reached = True

        if ~T_target_reached:  # if the target temperature has not been reached

            # if the rate is the same for next step, then the T ramp should also not change much
            # unless the continued T ramp will push the T significantly above T_target

            if rate_i[i] == rate_i[i+1]:
                T_final_next = T_i[i+1] + T_delta

                # if the projected T at the end of next step is below the target + threshold
                if T_final_next < (T_target+upper_thresh):
                    # keep ramp rate the same with some drift 'noise' to add variability
                    ramp_i.append(ramp_i[i]+getRandomDrift(rate_i[i+1]))

                else:
                    # change ramp so that it approaches T_target, with reasonable degC/min 'drift' rates
                    ramp_next = np.sign(
                        T_target-T_i[i+1])*abs(getRandomDrift(rate_i[i+1]))
                    ramp_i.append(ramp_next)

            # if the rate is not the same for next step, then a random T ramp should be selected
            # according to the previous T and current values

            else:
                ramp_next = getInitialRamp(T_i[i+1], rate_i[i+1])
                T_next_rand = T_i[i+1] + ramp_next*t_delta_i[i+1]

                # if the projected T at the end of next step is below the target + threshold
                if T_next_rand < (T_target+upper_thresh):
                    ramp_i.append(ramp_next)  # keep ramp rate the same

                # change ramp so that it approaches T_target, with reasonable degC/min 'drift' rates
                else:
                    ramp_next = np.sign(
                        T_target-T_i[i+1])*abs(getRandomDrift(rate_i[i+1]))
                    ramp_i.append(ramp_next)

        else:  # T target has already been reached, then temperature should drift around setpoint

            ramp_next = getRandomDrift(rate_i[i+1])
            ramp_i.append(ramp_next)

    # save parameters for jth set to lists
    list_ts.append(t_i)
    list_rates.append(rate_i)
    list_ramps.append(ramp_i)
    list_temps.append(T_i)
    list_SOCs.append(SOC_i)
    list_T_target.append(T_target)

#%% Create dataframe from parameters, hard-coding the variable names

df_times = 60*pd.DataFrame(data=np.array(list_ts),
                           columns=['t0', 't1', 't2', 't3', 't4', 't5', 't6', 't7', 't8']).drop(labels=['t0'], axis=1)

df_rates_all = pd.DataFrame(data=np.array(list_rates),
                            columns=['rate1', 'rate1__', 'rate2', 'rate2__', 'rate3', 'rate3__', 'rate4', 'rate4__', 'rate4___'])

df_rates = df_rates_all.copy()[['rate1', 'rate2', 'rate3', 'rate4']]

df_ramps = 1/60*pd.DataFrame(data=np.array(list_ramps),
                             columns=['ramp1', 'ramp2', 'ramp3', 'ramp4', 'ramp5', 'ramp6',
                                      'ramp7', 'ramp8', 'ramp9']).drop(labels=['ramp9'], axis=1)

df_temps = 273.15 + pd.DataFrame(data=np.array(list_temps),
                                 columns=['T0', 'T1', 'T2', 'T3', 'T4', 'T5', 'T6', 'T7', 'T8']).drop(labels=['T8'], axis=1)

df_target = 273.15 + pd.DataFrame(data=np.array(list_T_target), columns=['T_target'])

df_initialSOC = pd.DataFrame(data=np.array(
    list_SOCs).T[:][0], columns=['iSOC'])

df_all_params = pd.concat(
    [df_times, df_rates, df_ramps, df_temps, df_target, df_initialSOC], axis=1)

df_all_params = df_all_params.add_prefix('charge.')
#%% Functions for generating either correlated or random aging parameters

def getRandomAgingParams(n_simulations, bestValues, worstValues, names):
    """
    Goal: generate random aging parameter values for n_simulations assuming a 
    uniform random probability distribution for each parameter

    Inputs:
      n_simulations (float): number of sets of random values
      bestValues (list of floats): parameter values for no aging, beginning of life
      worstValues (list of floats): worst parameter values assumed possible in operational cell
      names (list of strings): parameter names to label DataFrame columns 

    Output:
      df_params (DataFrame): parameter values with shape (n_simulations,len(names))
    """
    randVals = np.random.rand(n_simulations, len(names))

    # matrix multiplication to duplicate the best and worst values by n_simulations
    bestVals = np.ones([n_simulations, 1])*np.array(bestValues).T
    worstVals = np.ones([n_simulations, 1])*np.array(worstValues).T

    # calculate parameter values
    paramValues = np.multiply(randVals, bestVals) + \
        np.multiply((1-randVals), worstVals)

    df_params = pd.DataFrame(data=paramValues, columns=names)

    return df_params

def getSOHbasedAgingParams(n_simulations, SOH_values, bestValues, worstValues, names):
    """
    Goal: generate aging parameter values correlated with the cell extent of aging.
    Here, SOH is used as a proxy for extent of aging, where 1 is the best parameters
    and 0.8 is end of life, corresponding to worst parameters

    Inputs:
      n_simulations (float): number of sets of random values
      SOH_values (list of floats): SOH values that specify the (mean, sd), for each set of parameters 
      bestValues (list of floats): parameter values for no aging, beginning of life
      worstValues (list of floats): worst parameter values assumed possible in operational cell
      names (list of strings): parameter names to label DataFrame columns 

    Output:
      df_params (DataFrame): parameter values with shape (n_simulations,len(names))

    Constants defined in 1st cell:
      SOH_min (float): SOH value that corresponds to 'worst' parameter values, ~0.8
      SOH_max (float): SOH value that corresponds to 'best' parameter values, 1
    """
    # get list of length n_simulations with approximately equal number of each SOH
    n_SOH = len(SOH_values)
    n_points_perSOH = n_simulations//n_SOH
    end_ind = np.ones(n_SOH)*n_points_perSOH
    end_ind[(n_SOH-1)] = end_ind[(n_SOH-1)]+n_simulations % n_SOH

    list_SOH = []
    for i in range(n_SOH):
        list_SOH = list_SOH+([SOH_values[i]]*int(end_ind[i]))

    # for list of SOH, return sampling of random nearby SOH, 1 for each aging parameter
    SOH_val = np.array(list_SOH).reshape([len(list_SOH), 1])
    n_params = len(names)

    # define standard deviation so that it increases with SOH (cell aging)
    # the size of this 2D array is n_simulationsxn_params
    sd = (1-SOH_val)/3*np.ones([1, n_params])
    mean = SOH_val*np.ones([1, n_params])
    rand_SOH_val = np.random.normal(mean, sd)

    # standardize SOH to values between 1 and 0, to allow scaling of parameters by best/worst param range
    scaleBy_0to1 = (rand_SOH_val-SOH_min)/(SOH_max-SOH_min)

    # matrix multiplication to duplicate the best and worst values by n_simulations
    bestVals = np.ones([n_simulations, 1])*np.array(bestValues).T
    worstVals = np.ones([n_simulations, 1])*np.array(worstValues).T

    # calculate parameter values by interpolating best and worst values with standardized SOH
    paramValues = np.multiply(scaleBy_0to1, bestVals) + \
        np.multiply((1-scaleBy_0to1), worstVals)

    df_params = pd.DataFrame(
        data=paramValues, columns=names)  # Create dataframe
    df_params['SOH'] = SOH_val

    return df_params

#%% Get aging parameters and export

if includeAgingParams:

    if agingCorrelatedSOH:
        df_aging_params = getSOHbasedAgingParams(n_conditions, SOH_values,
                                                 agingParamBestVal, agingParamWorstVal,
                                                 agingParamNames)
        fileName = '5June2023_params_SOHAging_finalParams_2000.csv'

    elif noAging:
        df_aging_params = getSOHbasedAgingParams(n_conditions, [1],
                                                 agingParamBestVal, agingParamBestVal,
                                                 agingParamNames)
        fileName = '15June2023_params_noAging_finalParams_1000.csv'
        
    else: #random aging parameters
        df_aging_params = getRandomAgingParams(n_conditions, agingParamBestVal,
                                               agingParamWorstVal, agingParamNames)
        # fileName = '19May2023_params_RandAging_correctedPorosity.csv'

    # add prefix to aging params to match COMSOL names
    df_aging_params = df_aging_params.add_prefix('aging.')

    # add df of random aging params to df_all_params
    df_all_withAging = pd.concat([df_all_params, df_aging_params], axis=1)

    # export df_all_params to csv
    df_all_withAging = df_all_withAging.round(decimals=7)

    # the anode Li loss can't exceed the available Li given the initial SOC,
    # so correct any of these cases
    an_LiLoss_series = df_all_withAging['aging.an_LiLoss']
    iSOC_series = df_all_withAging['charge.iSOC']
    notEnoughLi = an_LiLoss_series < iSOC_series
    lower_LiLoss_val = iSOC_series-0.02 * \
        np.random.rand(1)  # empirical correction
    an_LiLoss_modified = an_LiLoss_series.where(notEnoughLi, lower_LiLoss_val)

    # double check that correction is effective
    max(an_LiLoss_modified-iSOC_series)
    df_all_withAging['aging.an_LiLoss'] = an_LiLoss_modified

    # export file to CSV
    # df_all_withAging.to_csv(fileName)

else:
    df_all_params = df_all_params.round(decimals=4)
    # df_all_params.to_csv('15May2023_params_BOL_correctedPorosity.csv')


# double check that maximum temperature for all simulations doesn't exceed 70 degrees
# Occasionally, the algorithm can produce a non-physical high temperature value.
# In this case, generate a new set of parameters

maxTemp = 70
simMax_T = max(df_all_params.loc[:, 'charge.T0':'charge.T7'].max())-273.15
if simMax_T > maxTemp:
    print('WARNING, maximum temperature is: ' + str(np.round(simMax_T,3)) + ' degC')
    
else:
    print('Maximum temperature is: ' + str(np.round(simMax_T,3)) + ' degC')
    
#%% Visualize current and T profiles for some of the generated parameters 
fig, ax = plt.subplots(nrows=2, ncols=1, sharex=True)
fig.set_size_inches(8, 8)

for i in range(50, 200):
    # for i in range(n_conditions):
    plotSOC = list_SOCs[i]
    rates = list_rates[i]
    temps = list_temps[i]

    ax[0].step(plotSOC, rates, where='post')
    ax[1].plot(plotSOC[:], temps[:])

ax[0].set_xlim([0, 0.95])
ax[1].set_ylabel('T (C)')
ax[0].set_ylabel('C-rate')
ax[1].set_xlabel('SOC')