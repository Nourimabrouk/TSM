# -*- coding: utf-8 -*-
"""
Created on Wed Mar 18 22:53:48 2020

@author: sjoer
"""
import numpy as np
import pandas as pd
import math
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from scipy.stats import norm

def ParticleFilter(theta_tilde, w):
    
    N = len(theta_tilde.T) #number of draws
    T = len(theta_tilde)
    logp_yt_given_theta = np.zeros((T,N))
    logg_yt_given_theta = np.zeros((T,N))

    
    #NORMALIZE W
    w2 = np.zeros((T,N))
    for t in range(T):
        for n in range(N): 
            w2[t,n] = np.round(w[t,n]/ sum(w[t,:]),3)
    w = w2
    #draw new theta tilde: pick a theta_tilde(t,i) with weights w(t,i)
    theta_tilde_new = np.zeros((T,N))
    for t in range(T):
        vtheta                  = theta_tilde[t,:]
        vweights                = w[t,:]
        vweights[N-1]           = 1-sum(w[t,0:N-1]) #fix that weights sum up to 1
        for i in range(N): 
            theta_tilde_new[t,i]    = np.random.choice(vtheta, p=abs(vweights))
            
    ####Calculate new smoothed mean estimation    
            logp_yt_given_theta[:,i] = -0.5*(np.log(2*math.pi))- 0.5* theta_tilde_new[:,i] - 0.5*(y_tilde**2) * np.exp(- theta_tilde_new[:,i])
            logg_yt_given_theta[:,i] =   -0.5*(np.log(2*math.pi)) - 0.5 * np.log(A) - 0.5 *(y_star - theta_tilde_new[:,i])**2/A
#    
 
    m_sim = logp_yt_given_theta.sum(axis=0) - logg_yt_given_theta.sum(axis=0)
    #w_sim = np.exp(np.mean(m_sim))*np.exp(m_sim - np.mean(m_sim))    
    #x_hat = sum(theta_tilde * w_sim) / sum(w_sim)
 
    theta_hat_SPDK_new = np.zeros(T)
    for t in range(T):
        theta_hat_SPDK_new[t]    = sum(theta_tilde[t,:] * np.exp(m_sim - np.mean(m_sim)))/ sum(np.exp(m_sim - np.mean(m_sim)))

    ######PLOTS
    fig, axs = plt.subplots(1, 1, figsize=(15, 6), sharex=False, sharey=True)
    axs.plot(theta_tilde_new[1:T-1,])  
    axs.plot(H_bar[1:,], linewidth=3, color='black', label = 'Smoothed mean estimate')
    axs.plot(theta_hat_SPDK_new[1:,], linewidth=3, color='blue', label = 'Particle estimate')
    axs.legend()
    fig.suptitle("Particle filter for Gaussian distribution",fontsize = 20)
    plt.show()
          
    return theta_tilde_new, theta_hat_SPDK_new