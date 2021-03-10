import numpy as np
import pandas as pd
from scipy.optimize import minimize
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings("ignore")
import scipy.stats as stats
import random

def import_data(path):
    # import data and set date type to float32
    data = pd.read_csv(path,dtype='float32')
    # convert to np.array and flattens
    return pd.Series(data.to_numpy().flatten())


class KalmanFilter:

    is_initilized = False
    theta_hat = None

    def __init__(self, y, theta_ini, scaler = 100):

        # scaled back
        self.y = (y - np.mean(y))/scaler
        self.x = np.log(self.y**2)
        self.theta_ini = theta_ini
        


    def filtering(self, theta):
        
        x = self.x
        n = len(self.x)
        
        omega = theta[0]
        phi = theta[1]
        sigma2_eta = theta[2]

        # known sigma2_u_star
        sigma2_u_star = (np.pi**2) / 2
        
        
        # h eq. a
        h = np.zeros(n)
        # P 
        P = np.zeros(n)
        
        
        # initial values, unconditional mean and variance of h
        h[0] = omega/(1- phi)
        P[0] = sigma2_eta / (1 - phi**2)
        
        
        
        v = np.zeros(n)
        K = np.zeros(n)
        F = np.zeros(n)
        
        # filter
        for t in range(n-1):

            v[t] = x[t] - h[t] + 1.27 # prediction error in the Kalman filter
            F[t] = P[t] + sigma2_u_star #
            K[t] = P[t] / F[t] # Kalman gain
            h[t + 1] = omega + phi * h[t] + K[t] * v[t]

            # confidence bounds

            P[t + 1] = K[t] * sigma2_u_star + sigma2_eta

        
        # fill in the last value

        v[n-1] = x[n-1] - h[n-1] + 1.27
        F[n-1] = P[n-1] + sigma2_u_star
        K[n-1] = P[n-1] / F[n-1]    


        return v, F, h, P, K



        
    def llk_func(self, theta):
        '''
        calculate diffuse llk function
        '''
        n = len(self.x)  # number of observations

        v, F, _, _, _ = self.filtering(theta)

        llk = - (n / 2) * np.log(2 * np.pi) - (1 / 2) * \
            np.sum(np.log(F)) - (1 / 2) * np.sum(((v**2) / F))
        
        return -llk



    def fit(self, optim_method='SLSQP', display=False, max_iter = 300):
        '''
        minimize the negative llk, diffuse or concentrated
        llk_type: 'c' or 'd'
        '''

        print(f'Initial Theta value is set to {self.theta_ini}')
        optim_res = minimize(self.llk_func, self.theta_ini, method=optim_method, options={
                            'disp': display, 'maxiter': max_iter})
        
        self.theta_hat = optim_res.x

        return optim_res




    def particle_filtering(self, theta = None, seed = 1234):
        '''
        particle filter

        '''
        if theta is None:
            theta = self.theta_hat
            if theta is None:
                print('The parameters need to be estimated first.')

        # no stratified sampling, large N
        N = 10000
        omega = theta[0]
        phi = theta[1]
        sig2_eta = theta[2]

        xi = omega/(1-phi)

        y = self.y
        n = len(y)

        # set a random seed
        np.random.seed(seed)

        # alpha hat t|t
        a = np.zeros(n)

        # H tilde
        H = np.zeros((N,n))

        H[:,0] = np.random.normal(loc=0,scale = np.sqrt(sig2_eta/(1-phi**2)), size = N)

        for t in range(1,n):
            H[:,t] = np.random.normal(loc = phi*H[:,t-1],scale = np.sqrt(sig2_eta),size = N)
            #likeliehood
            w_tilde = stats.norm(0,np.sqrt(np.exp(xi)*np.exp(H[:,t]))).pdf(y[t])
            #Normalize
            w = w_tilde/np.sum(w_tilde)
            
            a[t] = np.sum(w * H[:,t])
            s = random.choices(range(N), weights = w, k = N)
            H = H[s,]

        # return alpha hat t|t
        return a


    def smoothing(self, theta = None):
        '''
        smoothing
        theta: use theta_hat
        '''

        # use the estimated parameters for smoothing
        if theta is None:
            theta = self.theta_hat
            if theta is None:
                print('The parameters need to be estimated first.')

        omega = theta[0]
        phi = theta[1]
        sigma2_eta = theta[2]
        
        v, F, h, P, K = self.filtering(theta)
        

        n = len(self.x)  # number of observations


        # initiation
        
        # for smoothed h
        h_hat = np.zeros(n)


        # no simulation atm
        V = np.zeros(n)

        # for smoothing
        L = np.zeros(n)
        r = np.zeros(n) # smoothing cumulant
        N = np.zeros(n) # smoothing varinace cumulant

        

        # n - 1, since len(0 ~ 99) = 100
        for t in range(n-1, -1, -1):

            L[t] = phi - K[t]

            # smoothing cumulant, r[n] = 0
            r[t-1] = v[t]/F[t] + L[t] * r[t]

            # smoothing varinace cumulant
            N[t-1] = 1/F[t] + (L[t]**2) * N[t]


            # alpha_hat
            h_hat[t] = h[t] + P[t] * r[t-1]
            V[t] = P[t] - (P[t]**2) * N[t-1]


        # confidence bounds
        upper_cb = h_hat + 1.645 * np.sqrt(V)
        lower_cb = h_hat - 1.645 * np.sqrt(V)


        return h_hat, upper_cb, lower_cb, v, F, h, P, K, L, r, N






if __name__ == '__main__':
    path= '../data/sv.dat'
    y = import_data(path)

    # initial values for optimization

    theta_ini = [0.001, # omega
                0.95, # phi
                0.5] # sigma2_eta


    kf = KalmanFilter(y, theta_ini)
    
    # estimating
    res = kf.fit()
    print(res)

    # # smoothing
    # h_hat, upper_cb, lower_cb, v, F, h, P, K, L, r, N = kf.smoothing()

    # fig = plt.figure(figsize=(20,5))
    # plt.scatter(x.index, x)
    # plt.plot(h_hat, color = 'black')
    # plt.show()

    a = kf.particle_filtering()
    fig = plt.figure(figsize=(20,5))  
    plt.plot(a, color = 'black')
    plt.show()
