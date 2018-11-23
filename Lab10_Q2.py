# -*- coding: utf-8 -*-
"""
Created on Wed Nov 21 09:50:16 2018

@author: brayd
"""

#PHY407 Lab10 Q2
#Monte-Carlo Integration to Estimate Volume of 10 dimensional hypersphere of unit radius
#From Newman Exercise 10.7

from random import uniform
import numpy as np

def f(r):
    '''
    A function that is the n-dimensional generalization of f(x,y) defined in Newman Ex. 10.7
    Input: r [array] of 1xD where D is the dimensionality of your hypersphere (in our case it'll be 10)
    Output: out [float] 1 if mod squared of r is less than or equal to 1, zero else
    '''
    normsquared = np.dot(r,r) #r*r = |r|^2 (n dimensional generalization of (x,y)*(x,y) = x^2+y^2)
    if normsquared <= 1:
        out = 1
    else:
        out = 0
    return out

D = 10 #dimensionality of our hypersphere
r = np.empty(D)

npoints = int(1e6) #a million points

I = 0
exp_f_squared = 0

#loop over the number of randomly sampled points we want
for i in range(npoints):
    #fill up r for each i
    for j in range(D):
        r[j] = uniform(-1,1) #random number between -1 and 1
    I += f(r)
    exp_f_squared += (f(r))**2

exp_f = I/npoints #<f>
exp_f_squared = exp_f_squared/npoints #<f^2>
#this is our approximation of the integral. V=2^D (volume of D-dimensional hypercube):
V = 2**D    
I = (V/npoints)*I

var_f = exp_f_squared - (exp_f)**2 #var f = <f^2> - <f>^2
std = V*np.sqrt(var_f)/np.sqrt(npoints) #standard deviation on integral according to eq. 10.32 Newman

print(\
'Volume of 10-dimension unit hypersphere estimated to be {0:.3f} with standard deviation {1:.2f}'.format(I, std))

