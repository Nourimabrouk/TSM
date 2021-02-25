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

def GetIniTheta():
    sigma_sq_eta    = 5
    phi             = 0.5
    omega           = 1  
   
    vTheta_ini      = [sigma_sq_eta, phi, omega]  
   
    sigma_sq_eta    = (0.001, 1000)
    phi             = (0.001, 0.999)
    omega           = (-1000, 1000)
    vTheta_bnds     = [sigma_sq_eta, phi, omega]
   
    return (vTheta_ini, vTheta_bnds)


def GetLogLikGauss(vTheta, data, sError):
    iT = len(data)        
    dLogLik = 0
    
    dSig, v, P, K, h     = KalmanFilterSV(data, vTheta)
    
    if(sError == 'Gauss'):
        for t in range(0,iT):
            logLikAtT = LogLikGauss(data[t], v[t] , dSig[t])
            dLogLik = dLogLik + logLikAtT
    return dLogLik

def LogLikGauss(dData, v, dSig):
    #if dSig > 0:
    return -(-0.5*math.log(2*math.pi) -0.5*math.log(dSig) - ((v)**2)/(2*dSig))


def KalmanFilterSV(data, vTheta):
    y       = data
    T       = len(data)

    sigma_sq_u      = (math.pi**2) / 2    
    mean_u          = -1.27
    
    sigma_sq_eta    = vTheta[0]
    phi             = vTheta[1]
    omega           = vTheta[2]
    
    # Kalman filtering
    T       = len(data)
    h       = np.zeros(T)
    P       = np.zeros(T)
    v       = np.zeros(T)
    F       = np.zeros(T)
    K       = np.zeros(T)
    
    h[0]    = omega / (1 - phi) #unconditional mean
    P[0]    = sigma_sq_eta / (1 - phi**2) # initial values 
    
    for t in range(T):
        v[t]    = y[t] -mean_u -h[t] 
        F[t]    = P[t] + sigma_sq_u
        K[t]    = phi*P[t]/F[t]
        #h_y[t]  = h[t] + P[t] * v[t] / F[t] 
        #P_y[t]  = P[t] - P[t]**2 / F[t]
        if t<T-1:
            h[t+1]  = omega + phi*h[t] + K[t]*v[t]
            P[t+1]  = (phi**2)*P[t] + sigma_sq_eta - (K[t]**2) * F[t]
    return F,v, P, K, h



def KalmanSmootherSV(data, vTheta):
    F,v, P, K, a = KalmanFilterSV(data, vTheta)
    
    phi             = vTheta[1]
    
    T           = len(v)
    alpha        = np.zeros(T)  #smoothed state
    N           = np.zeros(T)
    r           = np.zeros(T)
    V           = np.zeros(T) #smoothed state variance
    
    L           = phi - K
    N[T-1]      = 0
    r[T-1]      = 0
    
    for t in reversed(range(T)): #reversed loop
        if t>=0:
            N[t-1]      = (1/F[t])+L[t]**2 * N[t]
    
            if math.isnan(v[t]):
                r[t-1]      = r[t]
            else:                     
                r[t-1]      = (v[t]/F[t])+L[t]*r[t]
            
            
            V[t]                = P[t] - P[t]**2 * N[t-1]   
            alpha[t]            = a[t]+P[t]*r[t-1]
            

    return alpha

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

def mode_estimation(data, theta):
    sigma_sq_eta    = theta[0]
    phi             = theta[1]
    omega           = theta[2]  
    xi              = omega/ (1-phi)
    sigma           = np.exp(1/2* xi)
    
    y               = data[g_columnHeader]
    #y_tilde         = y- np.mean(y)/(sigma*np.exp(H_current))
    y_tilde         = y- np.mean(y)/(sigma)

    T               = len(y)
    sigma           = np.exp(1/2* xi)

    H_current         = np.ones(T)    
    A               = 2*np.ones(T)
    Z               = np.log(y_tilde**2)
    M               = 300

#    A[0]            = 2
#    Z[0]            = np.log()
    H               = np.zeros((T,M))
    for i in range(M): 
        #H[:,i]      = H_current
#        y_mean          = np.mean(y)
#        y_var           = sigma**2 * np.exp(H[t])
#        p               = -1/2 * np.log(2*math.pi) - 1/2 * np.log(sigma)-1/2*H[t] - 1/2 * y_tilde[t]**2 * np.exp(-H[t])
#        
        #y_tilde         = y- np.mean(y)/(sigma)
    

                
        P           = np.zeros(T)
        v           = np.zeros(T)
        F           = np.zeros(T)
        K           = np.zeros(T)
        
        H_current[0]    = omega / (1 - phi) #unconditional mean
        P[0]    = sigma_sq_eta / (1 - phi**2) # initial values         
        
        #Kalman Filter
        for t in range(T):
            v[t]    = Z[t] - H_current[t]
            F[t]    = P[t] + A[t]
            K[t]    = phi*P[t]/F[t]
            #h_y[t]  = h[t] + P[t] * v[t] / F[t] 
            #P_y[t]  = P[t] - P[t]**2 / F[t]
            if t<T-1:
                H_current[t+1]  = phi * H_current[t] + K[t]*v[t]
                P[t+1]  = phi**2 * P[t] + sigma_sq_eta - K[t]**2 * F[t]
        
       # print(v, F, K)
        #print(H_current)        
        
        #Smoothed Kalman Filter
        H_smooth    = np.zeros(T)
        P_smooth    = np.zeros(T)
        N           = np.zeros(T)
        r           = np.zeros(T)
        L           = phi - K
        N[T-1]      = 0
        r[T-1]      = 0
        for t in reversed(range(T)): #reversed loop
            if t>0:
                N[t-1]      = (1/F[t])+L[t]**2 * N[t]
        
                if math.isnan(v[t]):
                    r[t-1]      = r[t]
                else:                     
                    r[t-1]      = (v[t]/F[t])+L[t]*r[t]
                
                P_smooth[t]            = P[t] - P[t]**2 * N[t-1]   
                H_smooth[t]            = H_current[t]+P[t]*r[t-1]
                

                
        #print(H_smooth)
        H_current   = H_smooth
        
        A           = 2 * np.exp(H_current/y_tilde**2)
        Z           = H_current+1-1/2*A   
        
    
    theta_hat       = H_current
    A               = 2 * np.exp(theta_hat/y_tilde**2)
    y_star          = theta_hat +1-1/2*A  
    
    #check for convergence

        
    return theta_hat, A, sigma, T, sigma_sq_eta, phi, omega, H


def draw_theta(data, ML_theta):
    theta      = ML_theta
    M          = 10 #number of draws
    theta_hat, A, sigma, T, sigma_sq_eta, phi, omega, H  = mode_estimation(data, theta)
    sigma       = np.ones(T)*sigma
 
    
    theta_tilde = np.zeros((T,M)); theta_plus = np.zeros((T,M)); theta_plus_hat = np.zeros((T,M))
    eta_plus    = random_errors(0, sigma, M)
    u_plus      = random_errors(0, A, M)
    
    alpha_plus     = random_errors(np.log(sigma[0]), sigma, M) #initial value
    
    
    for i in range(M):

        for t in range(T): 
            if t>0:
                alpha_plus[t,i] = phi * alpha_plus[t-1,i] +eta_plus[t,i]
                theta_plus[t,i] = alpha_plus[t,i]  
                #theta_plus[t] = phi * theta_plus[t-1] +eta_plus
                y_plus      = theta_plus + u_plus
            

        H_current   = np.zeros(T)
        P           = np.zeros(T)
        v           = np.zeros(T)
        F           = np.zeros(T)
        K           = np.zeros(T)
        
        H_current[0]    = omega / (1 - phi) #unconditional mean
        P[0]    = sigma_sq_eta / (1 - phi**2) # initial values 
        
        for t in range(T):
            v[t]    = y_plus[t,i] - H_current[t]
            F[t]    = P[t] + A[t]
            K[t]    = phi*P[t]/F[t]
            #h_y[t]  = h[t] + P[t] * v[t] / F[t] 
            #P_y[t]  = P[t] - P[t]**2 / F[t]
            if t<T-1:
                H_current[t+1]  = phi * H_current[t] + K[t].T*v[t]
                P[t+1]  = phi**2 * P[t] + sigma_sq_eta - (K[t]**2).T * F[t]
                
        

        #Smoothed Kalman Filter
        H_smooth    = np.zeros(T)
        P_smooth    = np.zeros(T)
        N           = np.zeros(T)
        r           = np.zeros(T)
        L           = phi - K
        N[T-1]      = 0
        r[T-1]      = 0
        for t in reversed(range(T)): #reversed loop
            if t>0:
                N[t-1]      = (1/F[t])+L[t]**2 * N[t]
        
                if math.isnan(v[t]):
                    r[t-1]      = r[t]
                else:                     
                    r[t-1]      = (v[t]/F[t])+L[t].T*r[t]
                
                P_smooth[t]            = P[t] - P[t]**2 * N[t-1]   
                H_smooth[t]            = H_current[t]+P[t]*r[t-1]        
        
        theta_plus_hat[:,i] = H_smooth
        theta_tilde[:,i]    = theta_hat + theta_plus[:,i] - theta_plus_hat[:,i]

    plt.plot(theta_tilde)
    
#    ## construct densities
    ##Shapes kunnen niet met elkaar vermenigvuldigt worden
    
#    y                   = data[g_columnHeader]
#    logp_yt_given_omega = -0.5*(np.log(2*math.pi*sigma) + theta_tilde.T + np.exp(- theta_tilde) * (y**2/sigma))
#    logg_yt_given_omega = -1/2*np.log(2*math.pi*A) + theta_tilde**2 / A
#    
#    m_sim = sum(logp_yt_given_omega) - sum(logg_yt_given_omega)
#    w_sim = np.exp(np.mean(m_sim))*np.exp(m_sim - np.mean(m_sim))
#    
#    x_hat = sum(theta_tilde * w_sim) / sum(w_sim)
    
    
    return theta_tilde, H #, plt.plot(theta_tilde)



def random_errors(mu, var, M):
    T= len(var)
    error = np.zeros((T,M))
    for t in range(T):
        for m in range(M):
            error[t,m]   = np. random. normal(mu,var[t])
        
    return error

def LoadData():
    global g_columnHeader, g_xColumnHeader
    global g_x, g_sv
    g_columnHeader = 'Financial returns'
    g_xColumnHeader = 'Transformed financial returns'
    data =  np.genfromtxt('sv.dat', skip_header=1)
    g_sv = pd.DataFrame(data)/100
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
    global g_sv, g_x, vTheta_ML
    
    LoadData()
    
    vTheta_ML, llik = GetOptKalman("Gauss")
    alpha = KalmanSmootherSV(g_x[g_xColumnHeader], vTheta_ML)
    theta_tilde, H = draw_theta(g_sv, vTheta_ML)


if __name__ == "__main__":
    main()