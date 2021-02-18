# -*- coding: utf-8 -*-
"""
Code for assignment two of the Time Series Models course.
Created on 25-02-2020 by Job
"""

import numpy as np
import pandas as pd
import math
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import plotting_ass_2
from scipy.stats import norm
from scipy.optimize import fmin_slsqp

##############################################################################
####KALMAN SECTION
def InitialValues(model):
    if model == "KalmanFilter": 
        a_ini   = 0
        P_ini   = 10**7 
        theta   = [a_ini, P_ini]
    return theta

def KalmanFilter(data, sigma_sq_eps2, sigma_sq_eta):
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
        F[t]    = P[t] + sigma_sq_eps2
        
        if math.isnan(y[t]):
            K[t]    = 0
            a_y[t]  = a[t]
            P_y[t]  = P[t]
            a[t+1]  = a[t] 
            P[t+1]  = P[t] + sigma_sq_eta
        else:
            K[t]    = P[t]/F[t]
            a_y[t]  = a[t] + K[t] * v[t]
            P_y[t]  = P[t]*(1 - K[t])
            
            if t< T-1:
                a[t+1]  = a[t] + K[t]*v[t]
                P[t+1]  = P[t]*(1 - K[t]) + sigma_sq_eta            
                
    a_lb    = a - 1.645*np.sqrt(P) #upper and lower bounds of a    
    a_ub    = a + 1.645*np.sqrt(P)
    
    DF_Kalman           = pd.DataFrame([a, P, v, F, K, a_y, P_y, a_lb, a_ub]).T
    DF_Kalman.columns   = ["a", "P", "v", "F", "K", "a_y", "P_y", "a_lb", "a_ub"]
    return DF_Kalman

def GetOptKalman(sError):
    IniTheta = GetIniTheta()
    vIniTheta = IniTheta[0]
    vIniThetaBounds = IniTheta[1]
    results     = fmin_slsqp(GetLogLikGauss, vIniTheta, args = (g_x[g_xColumnHeader],sError), bounds = vIniThetaBounds, full_output= True)
    vTheta_ML, llik   = results[0], -results[1]
    print("\nTheta estimate:", vTheta_ML)
    print("\nLog likelihood value:", llik)
    return vTheta_ML, llik
GetOptKalman("Gauss")




def GetIniTheta():
    sigma_sq_eta    = 5
    phi             = 0.5
    omega           = 1  
   
    vTheta_ini      = [sigma_sq_eta, phi, omega]  
   
    sigma_sq_eta    = (0.001, 1000)
    phi             = (0, 1)
    omega           = (-1000, 1000)
    vTheta_bnds     = [sigma_sq_eta, phi, omega]
   
    return (vTheta_ini, vTheta_bnds)


def GetLogLikGauss(vTheta, data, sError):
    iT = len(data)        
    dLogLik = 0
    
    dSig, mu     = KalmanFilterSV(data, vTheta)
    
    if(sError == 'Gauss'):
        for t in range(1,iT-1):
            logLikAtT = LogLikGauss(data[t], mu[t] , dSig[t])
            dLogLik = dLogLik + logLikAtT
    return dLogLik

def LogLikGauss(dData, mu, dSig):
    if dSig > 0:
        return -(- 0.5*math.log(2*math.pi) - 0.5*math.log(dSig) - ((dData-mu)**2)/(2*dSig))

def KalmanFilterSV(data, vTheta):
    y       = data
    
    sigma_sq_eps       = (math.pi) **2 / 2    
    
    sigma_sq_eta = vTheta[0]
    phi = vTheta[1]
    omega = vTheta[2]
    mean_u = -1.27
    
    # Kalman filtering
    T       = len(data)
    h       = np.zeros(T)
    P       = np.zeros(T)
    v       = np.zeros(T)
    F       = np.zeros(T)
    K       = np.zeros(T)
    
    h_y     = np.zeros(T) #conditional on y
    P_y     = np.zeros(T) 
    
    h[0]    = 0
    P[0]    = sigma_sq_eta # initial values 
    
    for t in range(T):
        v[t]    = y[t] - h[t] - mean_u
        F[t]    = P[t] + sigma_sq_eps
        K[t]    = phi*P[t]/F[t]
        h_y[t]  = h[t] + P[t] * v[t] / F[t]
        P_y[t]  = P[t] - P[t]**2 / F[t]
        if t<T-2:
            h[t+1]  = phi * h_y[t] + omega
            P[t+1]  = phi**2 * P_y[t] + sigma_sq_eta
    return P,h


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
    
    return alpha, N, r, V, alpha_lb, alpha_ub

def LoadData():
    global g_columnHeader, g_xColumnHeader
    global g_x, g_sv
    g_columnHeader = 'Financial returns'
    g_xColumnHeader = 'Transformed financial returns'
    data =  np.genfromtxt('sv.dat', skip_header=1)
    g_sv = pd.DataFrame(data)
    g_sv.columns = [g_columnHeader]
    svStatistics = g_sv.describe()
    print("\nSome basic statistics on the financial returns data:\n",svStatistics)
    svMean = svStatistics[g_columnHeader][1]
    plotting_ass_2.PlotReturns(g_sv,g_columnHeader)
    g_x = g_sv.copy()
    g_x.columns = [g_xColumnHeader]
    g_x[g_xColumnHeader] = np.log((g_x[g_xColumnHeader] - svMean)**2)
    xStatistics = g_x.describe()
    print("\nSome basic statistics on the transformed data:\n",xStatistics)
    plotting_ass_2.PlotReturns(g_x,g_xColumnHeader)
    return

##############################################################################
####Main Function
def main():
#Define variables for global environment    
    global g_sv, g_x
    
    LoadData()
    
    vTheta_ML, llik = GetOptKalman("Gauss")
    
if __name__ == "__main__":
    main()