#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 18 16:10:43 2020

@author: jobdenotter
"""
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

def PlotReturns(svData, label):
    print("\nThe following figure shows a plot of the "+label)
    plt.figure(figsize=(15,10))
    plt.plot(svData, label = label)
    plt.legend()
    return plt.show()

def Plot2Lines(kalman, kalmanLabel,data, dataLabel):
    print("\nThe following plot shows the original data "+kalmanLabel +" with filtered kalman state a")
    plt.figure(figsize=(15,10))
    plt.plot(kalman, color = "black", label = kalmanLabel)
    plt.plot(data, color = "lightblue", label = dataLabel)
    plt.legend()
    return plt.show()

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

#def PlotKalmanSV():
#    print("\n ")
#    plt.figure(figsize=(20,10))
#    plt.plot(h_smoothing, color = "black", label = kalmanLabel)
#    plt.plot(data, color = "lightblue", label = dataLabel)
#    plt.legend()
#    return plt.show()