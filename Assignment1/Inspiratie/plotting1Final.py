#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 18 16:10:43 2020

@author: jobdenotter
"""
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
import seaborn as sns
from statsmodels.tsa.stattools import acf



def PlotForecast(DF_forecast, DF_Kalman, NileData_array, yearIncludingForecast, T):
    a = np.array(DF_Kalman["a"])[1:]
    P = np.array(DF_Kalman["P"])[1:]
    F = np.array(DF_Kalman["F"])[1:]
    
    forecast_a = np.array(DF_forecast["a"])
    forecast_a_lb = np.array(DF_forecast["a_lb"])
    forecast_a_ub = np.array(DF_forecast["a_ub"])
    forecast_P = np.array(DF_forecast["P"])
    forecast_F = np.array(DF_forecast["F"])
    
    total_a = np.concatenate((a,forecast_a))
    total_P = np.concatenate((P,forecast_P))
    total_F = np.concatenate((F,forecast_F))
    
    fig, axs = plt.subplots(2, 2, figsize=(15, 15), sharex=False, sharey=False)
    axs[0,0].plot(yearIncludingForecast[1:], total_a, color='black')
    axs[0,0].plot(yearIncludingForecast[T:], forecast_a_lb, color='grey')
    axs[0,0].plot(yearIncludingForecast[T:], forecast_a_ub, color='grey') 
    axs[0,0].plot(yearIncludingForecast[:T],NileData_array, 'o', color='black')
    
    axs[0,1].plot(yearIncludingForecast[1:], total_P, color='black')
    
    axs[1,0].plot(yearIncludingForecast[1:], total_a, color='grey')
    
    axs[1,1].plot(yearIncludingForecast[1:], total_F, color='black')
    fig.suptitle("Forecasting",fontsize = 20)
    return plt.show()


def PlotKalman(NileData, DF_Kalman):
    fig, axs = plt.subplots(2, 2, figsize=(15, 15), sharex=False, sharey=False)   
    
    axs[0,0].plot(NileData, 'o', color='grey')
    axs[0,0].plot(DF_Kalman.a[1:], color='black')
    axs[0,0].plot(DF_Kalman.a_lb[1:], color='red')
    axs[0,0].plot(DF_Kalman.a_ub[1:], color='red')    
    
    axs[0,1].plot(DF_Kalman.P[1:], color='black')
    
    axs[1,0].plot(DF_Kalman.v[1:], color='black')
    axs[1,0].axhline(y=0.0, color='r', linestyle='-')
    axs[1,1].plot(DF_Kalman.F[1:], color='black')
    
    axs[0,0].set_ylim([450, 1400])
    axs[0,1].set_ylim([5000, 17500])
    axs[1,0].set_ylim([-450, 450])
    axs[1,1].set_ylim([20000, 32500])
    
    fig.suptitle("Kalman Filter",fontsize = 20)

    return plt.show()

def PlotSmoothedState(NileData,DF_SmoothedState):
    fig, axs = plt.subplots(2, 2, figsize=(15, 15), sharex=False, sharey=False)   
    
    axs[0,0].plot(NileData, 'o', color='grey')
    axs[0,0].plot(DF_SmoothedState.alpha[1:-1], color='black')
    axs[0,0].plot(DF_SmoothedState.alpha_lb[1:-1], color='red')
    axs[0,0].plot(DF_SmoothedState.alpha_ub[1:-1], color='red')    
    axs[0,1].plot(DF_SmoothedState.V[1:], color='black')
    
    axs[1,0].plot(DF_SmoothedState.r[:-1], color='black')
    axs[1,0].axhline(y=0.0, color='r', linestyle='-')
    axs[1,1].plot(DF_SmoothedState.N[:-1], color='black')

    
    axs[0,0].set_ylim([450, 1400])
    axs[0,1].set_ylim([2300, 4100])
    axs[1,0].set_ylim([-0.036, 0.024])
    
    fig.suptitle("Smoothed State",fontsize = 20)

    return plt.show()

def PlotMissingObs(NileData_MissingObs, DF_MissingObs_Kalman, DF_MissingObs_SmoothedState):

    fig, axs = plt.subplots(2, 2, figsize=(15, 15), sharex=False, sharey=False)
    
    axs[0,0].plot(DF_MissingObs_Kalman.a[1:], color='black')
    axs[0,0].plot(NileData_MissingObs[1:], color='grey')
        
    axs[0,1].plot(DF_MissingObs_Kalman.P[1:], color='black')
    
    axs[1,0].plot(NileData_MissingObs[1:], color='grey')
    axs[1,0].plot(DF_MissingObs_SmoothedState.alpha[1:], color='black')

    axs[1,1].plot(DF_MissingObs_SmoothedState.V[1:], color='black')
    
    fig.suptitle("Filtering and Smoothing output when observations are missing",fontsize = 20)    
    return plt.show()

def PlotDistSmo(DF_DisturbanceSmoothing):
    
    eps = np.array(DF_DisturbanceSmoothing["eps"])
    sd_eps = np.array(DF_DisturbanceSmoothing["sd_eps"])
    eta = np.array(DF_DisturbanceSmoothing["eta"])
    sd_eta = np.array(DF_DisturbanceSmoothing["sd_eta"])
    
    fig, axs = plt.subplots(2, 2, figsize=(15, 15), sharex=False, sharey=False)
    
    axs[0,0].plot(eps, color='black')
    axs[0,1].plot(sd_eps, color='black')
    axs[1,0].plot(eta, color='black')
    axs[1,1].plot(sd_eta, color='black')
    
    axs[0,0].axhline(y=0.0, color='r', linestyle='-')
    axs[1,0].axhline(y=0.0, color='r', linestyle='-')
    
    fig.suptitle("Output of disturbance smoothing recursion",fontsize = 20)
    
    return plt.show()

def PlotStandardPredictionErrors(DF_st_err):
    fig, axs = plt.subplots(2, 2, figsize=(15, 15), sharex=False, sharey=False)   
    
    axs[0,0].plot(DF_st_err['stv'], '-', color='k')
    axs[0,1].hist(DF_st_err['stv'], 14, density=1, color='w', edgecolor='k')
    
    sns.distplot(DF_st_err['stv'], hist=False, rug=False, color = "black", ax = axs[0,1])
    axs[0,1].set_ylabel('')    
    axs[0,1].set_xlabel('')
    

    stats.probplot(DF_st_err['stv'], dist="norm", plot=axs[1,0])
    axs[1,0].get_lines()[0].set_linestyle('-')
    axs[1,0].get_lines()[0].set_marker(None)
    axs[1,0].get_lines()[0].set_color('k')
    axs[1,0].get_lines()[1].set_color('k')
    axs[1,0].set_ylabel('')    
    axs[1,0].set_xlabel('')
    axs[1,0].set_title('')
    
    nlags = 10
    corr = acf(DF_st_err['stv'], nlags=nlags)
    axs[1,1].bar(x=np.arange(1,nlags+1), height=corr[1:], color='grey', edgecolor='k')
    axs[1,1].set_ylim([-1,1])

    


    #fig.suptitle("Standardized residuals",fontsize = 20)

    return plt.show()

def PlotSmoothResiduals(DF_st_smooth_res):
    fig, axs = plt.subplots(2, 2, figsize=(15, 15), sharex=False, sharey=False)   
    
    axs[0,0].plot(DF_st_smooth_res['u_star'], '-', color='k')
    axs[0,1].hist(DF_st_smooth_res['u_star'], 14, density=1, color='w', edgecolor='k')
    sns.distplot(DF_st_smooth_res['u_star'], hist=False, rug=False, ax = axs[0,1], color = 'k')
    axs[0,1].set_ylabel('')    
    axs[0,1].set_xlabel('')
    
    axs[1,0].plot(DF_st_smooth_res['r_star'], '-', color='k')
    axs[1,1].hist(DF_st_smooth_res['r_star'], 14, density=1, color='w', edgecolor='k')
    sns.distplot(DF_st_smooth_res['r_star'], hist=False, rug=False, ax = axs[1,1], color = 'k')
    axs[1,1].set_ylabel('')    
    axs[1,1].set_xlabel('')

    #fig.suptitle("Standardized smoothed residuals",fontsize = 20)
    return plt.show()