# -*- coding: utf-8 -*-
"""
Created on Wed Nov 21 10:24:28 2018

@author: brayd
"""

from random import random
import numpy as np
import matplotlib.pyplot as plt


#part a: 
def integrand(x):
    return x**(-1/2)/(1+np.exp(x))

N = 10000 #number of sampled points
npasses = 100 #number of times we will estimate integral with each method
#==============================================================================
#regular mean value method

Imv = np.zeros(npasses) #array to hold integral estimations for mean value method

#loop over number of estimations
for t in range(npasses):
    #perform mean value monte-carlo integration
    for i in range(N):
        xval = random() #random value in [0,1)
        Imv[t] += integrand(xval) #add value of integrand at sampled point to the sum
Imv = (1/N)*Imv #add coefficient of (b-a)/N, here b=1, a=0

#plot histogram
plt.figure()
plt.title('Mean Value estimation distribution for $\int_0^1 x^{-1/2}/(1+exp(x))dx$')
plt.xlabel('Estimation value')
plt.ylabel('Occurance')
plt.hist(Imv, 10,range=[0.8, 0.88])
plt.show()
#==============================================================================
#importance sampling method

Iis = np.zeros(npasses) #array to hold integral estimations for importance sampling

def p(x):
    '''
    Distribution we wish to sample from
    '''
    return 1/(2*x**(1/2))

def transform(z):
    x = z**2
    return x

#loop over number of estimations
for t in range(npasses):
    #perform importance sampling monte-carlo integration
    for i in range(N):
        z = random() #random value in [0,1) from uniform distribution
        xval = transform(z) #transform to sample x from p(x) distribution 
        #add value of integrand at sampled point to the sum according to Newman eq. 10.42
        Iis[t] += integrand(xval)/p(xval)
Iis = (1/N)*Iis #add coefficient of (b-a)/N, here b=1, a=0

#plot histogram
plt.figure()
plt.title('Importance sampling estimation distribution for $\int_0^1 x^{-1/2}/(1+exp(x))dx$')
plt.xlabel('Estimation value')
plt.ylabel('Occurance')
plt.hist(Iis, 10,range=[0.8, 0.88])
plt.show()
#==============================================================================
#==============================================================================
#part b:
def integrand_b(x):
    return np.exp(-2*abs(x-5))

#we'll reuse npasses and N from a
    
#==============================================================================
#regular mean value method

Imv_b = np.zeros(npasses) #array to hold integral estimations for mean value method

#loop over number of estimations
for t in range(npasses):
    #perform mean value monte-carlo integration
    for i in range(N):
        xval = random() #random value in [0,1)
        Imv_b[t] += integrand_b(xval) #add value of integrand at sampled point to the sum
Imv_b = (1/N)*Imv_b #add coefficient of (b-a)/N, here b=1, a=0

#plot histogram
plt.figure()
plt.title('Mean Value estimation distribution for $\int_0^{10} exp(-2|x-5|)dx$')
plt.xlabel('Estimation value')
plt.ylabel('Occurance')
plt.hist(Imv_b, 10)
plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
plt.show()

#==============================================================================
#importance sampling method

Iis_b = np.zeros(npasses) #array to hold integral estimations for importance sampling

def p_b(x):
    '''
    Distribution we wish to sample from (Gaussian)
    '''
    return 1/(np.sqrt(2*np.pi))*np.exp(-(x-5)**2/2)
#loop over number of estimations
for t in range(npasses):
    #perform importance sampling monte-carlo integration
    for i in range(N):
        xval = np.random.normal(loc = 5.0, scale = 1.0)  #random value from p_b(x Gaussian distribution
        #add value of integrand at sampled point to the sum according to Newman eq. 10.42
        Iis_b[t] += integrand_b(xval)/p_b(xval)
Iis_b = (1/N)*Iis_b #add coefficient of (b-a)/N, here b=1, a=0

#plot histogram
plt.figure()
plt.title('Importance sampling estimation distribution for $\int_0^{10} exp(-2|x-5|)dx$')
plt.xlabel('Estimation value')
plt.ylabel('Occurance')
plt.hist(Iis_b, 10)
plt.show()
