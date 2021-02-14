# -*- coding: utf-8 -*-
"""
Created on Thu Feb 13 12:16:48 2020

@author: sjoer
"""

import numpy as np
import pandas as pd
import math
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from scipy.stats import norm

from plotting import PlotForecast, PlotKalman, PlotSmoothedState, PlotMissingObs


##############################################################################
####KALMAN SECTION
#Given parameters
#below values are variances!
sigma_eps       = 15099
sigma_eta       = 1469.1


def InitialValues(model):
    if model == "KalmanFilter": 
        a_ini   = 0
        P_ini   = 10**7 
        theta   = [a_ini, P_ini]
        
    return theta


def KalmanFilter(data, sigma_eps, sigma_eta):
    name    = "KalmanFilter"
    y       = data
    T       = len(data)
    a       = np.zeros(T)
    P       = np.zeros(T)
    v       = np.zeros(T)
    F       = np.zeros(T)
    K       = np.zeros(T)
    
    a_y     = np.zeros(T) #conditional on y
    P_y     = np.zeros(T) 
    
    a[0], P[0]    = InitialValues(name) # intitial values 
    
    for t in range(T):
        v[t]    = y[t] - a[t]
        F[t]    = P[t] + sigma_eps
        
        if math.isnan(y[t]):
            K[t]    = 0
            a_y[t]  = a[t]
            P_y[t]  = P[t]
            a[t+1]  = a[t] 
            P[t+1]  = P[t] + sigma_eta
        else:
            K[t]    = P[t]/F[t]
            a_y[t]  = a[t] + K[t] * v[t]
            P_y[t]  = P[t]*(1 - K[t])
            
            if t< T-1:
                a[t+1]  = a[t] + K[t]*v[t]
                P[t+1]  = P[t]*(1 - K[t]) + sigma_eta            
                

            
    a_lb    = a - 1.645*np.sqrt(P) #upper and lower bounds of a    
    a_ub    = a + 1.645*np.sqrt(P)
    
    DF_Kalman           = pd.DataFrame([a, P, v, F, K, a_y, P_y, a_lb, a_ub]).T
    DF_Kalman.columns   = ["a", "P", "v", "F", "K", "a_y", "P_y", "a_lb", "a_ub"]
    DF_Kalman.index     = np.arange(1871, 1971); DF_Kalman.index.name = 'Year'
    return DF_Kalman


def SmoothedState(a, P, v, F, K):

    T           = len(v)
    alpha       = np.zeros(T)  
    N           = np.zeros(T)
    r           = np.zeros(T)
    V           = np.zeros(T) #smoothed state variance
    
    L           = 1 - K
    N[T-1]      = 0
    r[T-1]      = 0
    
    for t in reversed(range(T)): #reversed loop
        if t>0:
            N[t-1]      = (1/F[t])+L[t]**2 * N[t]
    
            if math.isnan(v[t]):
                r[t-1]      = r[t]
            else:                     
                r[t-1]      = (v[t]/F[t])+L[t]*r[t]
            
            
            V[t]                = P[t] - P[t]**2 * N[t-1]   
            alpha[t]            = a[t]+P[t]*r[t-1]
            

    alpha_lb    = alpha - 1.645*np.sqrt(V) #upper and lower bounds of a    
    alpha_ub    = alpha + 1.645*np.sqrt(V)

    DF_SmoothedState            = pd.DataFrame([alpha, N, r, V, alpha_lb, alpha_ub]).T
    DF_SmoothedState.columns    = ["alpha", "N", "r", "V", "alpha_lb", "alpha_ub"]
    DF_SmoothedState.index      = np.arange(1871, 1971)
    DF_SmoothedState.index.name = 'Year' 
    
    return DF_SmoothedState

def DisturbanceSmoothing(F,v,K,r,N):
#calculations
    u = 1/F * v - K * r
    eps = (sigma_eps) * u 
    D = 1/F + np.square(K) * N
    var_eps = (sigma_eps) - np.square(sigma_eps) * D
    eta = sigma_eta * r
    var_eta = sigma_eta - np.square(sigma_eta) * N
    sd_eta = np.sqrt(var_eta)
    sd_eps = np.sqrt(var_eps)
    year = np.arange(1871, 1971)
    
#plots
    fig, axs = plt.subplots(2, 2, figsize=(15, 15), sharex=False, sharey=False)
    
    axs[0,0].plot(year, eps, color='black')
    axs[0,1].plot(year, sd_eps, color='black')
    axs[1,0].plot(year, eta, color='black')
    axs[1,1].plot(year, sd_eta, color='black')
    
    axs[0,0].axhline(y=0.0, color='r', linestyle='-')
    axs[1,0].axhline(y=0.0, color='r', linestyle='-')
    
    fig.suptitle("Output of disturbance smoothing recursion",fontsize = 20)
    
    return plt.show()

def Forecasting(numberOfForecastSteps,startYear,T):
    yearIncludingForecast = np.arange(startYear, startYear + T + numberOfForecastSteps)

    forecast_a = np.zeros(numberOfForecastSteps)
    forecast_P = np.zeros(numberOfForecastSteps)
    forecast_F = np.zeros(numberOfForecastSteps)

    forecast_a[0] = np.array(DF_Kalman["a"])[-1]
    forecast_P[0] = np.array(DF_Kalman["P"])[-1] + sigma_eta
    forecast_F[0] = forecast_P[0] + sigma_eps
    for j in (range(numberOfForecastSteps-1)): #reversed loop
        forecast_a[j+1] = forecast_a[j]
        forecast_P[j+1] = forecast_P[j] + sigma_eta
        forecast_F[j+1] = forecast_P[j+1] + sigma_eps
        
    forecast_a_lb    = forecast_a - 0.675*np.sqrt(forecast_F)   
    forecast_a_ub    = forecast_a + 0.675*np.sqrt(forecast_F)

    DF_forecast = pd.DataFrame([forecast_a, forecast_P, forecast_F,forecast_a_lb,forecast_a_ub]).T
    DF_forecast.columns = ["a", "P", "F", "a_lb", "a_ub"]
    DF_forecast.index = np.arange(1971, 1971 + numberOfForecastSteps)
    DF_forecast.index.name = 'Year'
    
    PlotForecast(DF_forecast, DF_Kalman, NileData_array, yearIncludingForecast, T)
    return DF_forecast
    
####Data Preperation

def LoadData():
    data =  np.genfromtxt('Nile.dat', skip_header=1)
    NileData_array = data
    NileData = pd.DataFrame(data)
    NileData.columns = ['Annual Flow Nile']
    NileData.index = np.arange(1871, 1971)
    NileData.index.name = 'Year'
    
    return NileData, NileData_array 
    
def LoadDataMissing():    ####Missing observations
    data =  np.genfromtxt('Nile.dat', skip_header=1)
    NileData_array2 = data

    
    NileData_MissingObs_array           = NileData_array2
    NileData_MissingObs_array[20:40]    = math.nan
    NileData_MissingObs_array[60:80]    = math.nan
    NileData_MissingObs                 = pd.DataFrame(NileData_MissingObs_array)
    NileData_MissingObs.index           = np.arange(1871, 1971); NileData_MissingObs.index.name = 'Year'

    
    return NileData_MissingObs, NileData_MissingObs_array


##############################################################################
####Main Function
def main():

#Define variables for global environment    
    global NileData, NileData_array, NileData_MissingObs, NileData_MissingObs_array, DF_Kalman, DF_SmoothedState, DF_MissingObs_Kalman, DF_MissingObs_SmoothedState
    #Load Nile Data and Missing Observations
    NileData, NileData_array = LoadData()
    NileData_MissingObs, NileData_MissingObs_array  = LoadDataMissing()
    
    ####KALMAN FILTER
    DF_Kalman = KalmanFilter(NileData_array, sigma_eps, sigma_eta)
    PlotKalman(NileData, DF_Kalman)   
    
    ####SMOOTHED STATE
    #Input
    a, P, v, F, K   = DF_Kalman['a'].values, DF_Kalman['P'].values, DF_Kalman['v'].values,  DF_Kalman['F'].values, DF_Kalman['K'].values   
    #Smoothed State Output
    DF_SmoothedState    = SmoothedState(a,P,v,F,K)
    PlotSmoothedState(NileData,DF_SmoothedState)
    
    ####DISTURBANCE SMOOTHING
    #Input
    F,v,K,r,N = DF_Kalman['F'].values, DF_Kalman['v'].values, DF_Kalman['K'].values, DF_SmoothedState['r'].values, DF_SmoothedState['N'].values
    DisturbanceSmoothing(F,v,K,r,N)

    ####MISSING OBSERVATIONS
    #Missing observation DataFrame for Kalman & SmoothedState:
    DF_MissingObs_Kalman                = KalmanFilter(NileData_MissingObs_array, sigma_eps, sigma_eta)
    DF_MissingObs_SmoothedState         = SmoothedState(DF_MissingObs_Kalman['a'].values, DF_MissingObs_Kalman['P_y'].values, DF_MissingObs_Kalman['v'].values, DF_MissingObs_Kalman['F'].values, DF_MissingObs_Kalman['K'].values)
    PlotMissingObs(NileData_MissingObs, DF_MissingObs_Kalman, DF_MissingObs_SmoothedState)

    ####FORECASTING
    Forecasting(30, 1871,100)
    
if __name__ == "__main__":
    main()