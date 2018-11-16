# -*- coding: utf-8 -*-
"""
Created on Wed Nov 14 09:46:30 2018

@author: brayd
"""

import numpy as np
from scipy.constants import m_e, hbar
#instantiate some values of constants 
L = 1e-8 #length of box in meters
N = 1000 #spatial slices
a = L/N #size of grid boxes
M = m_e #mass = mass of electron in kg

#instantiate emtpy vector for psi(x,0)
psi = np.zeros(N+1, dtype = complex)
#define the function that psi(x,0) is represented by to later fill up the vector
def fun_psi_0(x):
    '''
    Calculate the wavefunction at t=0, for a value of x. 
    x should be a real valued function between zero and L for this to make
    sense, physically. 
    '''
    x_o = L/2 #mean of Guassian part of the wavefunction in meters
    sigma = 1e-10 #stdev of Gaussian part of the wavefunction in meters
    k = 5e10 #m^{-1} 
    out = np.exp(-(x - x_o)**2/(2*sigma**2))*np.exp(1j*k*x)
    #NB: we forget about a normalization constant because it doesn't affect the solution
    return out
#fill up the array, leaving the boundary conidtions in place
x = np.linspace(0, L, N+1)
psi[:] = fun_psi_0(x)
psi[0] = 0
psi[N] = 0
#normalize
norm = np.trapz(abs(psi)**2)
psi = psi/norm
    
#define v := B*psi(x,0), B as given in Newman pg.440
h = 1e-18 #time step in seconds
b1 = 1 - h*1j*hbar/(2*M*a**2)
b2 = h*1j*hbar/(4*M*a**2)
v = np.zeros(N, dtype = complex)

#calculate v = B.psi(x,0)
v = b1*psi[1:N] + b2*(psi[2:N+1] + psi[0:N-1])
    
#now to solve the system Ax = v
#x will be time evolved value of psi, A is as given in Newman pg.440
from banded import banded
a1 = 1 + h*1j*hbar/(2*M*a**2)
a2 = -h*1j*hbar/(4*M*a**2)
A = np.empty([3,N], dtype = complex)
A[0,:] = a2
A[1,:] = a1
A[2,:] = a2
psi_t = np.zeros(N+1, dtype = complex) #psi at later time t
#Calculate 1 time step: 
psi[1:N] = banded(A, v, 1, 1) #solve for x using Gauss elim and backsub for tridiagnonal matrix
#the above psi_t is the wavefunction on the domain x in [0,L] for t=h (h time step)

#Extend program to perform repeated steps and make an animation 

steps = 2000
from pylab import figure, plot, xlabel, ylabel, title, clf, pause, draw

figure()
plot(x, psi.real)
xlabel('$x$ m')
ylabel('$\psi(x,t)$')
title('Wavefunction at time t = 0')

#instantiate emtpy vector for psi(x,t)
psi_t = np.zeros(N+1, dtype = complex)
#fill up the array, leaving the boundary conidtions in place
psi_t[:] = fun_psi_0(x)
psi_t[0] = 0
psi_t[N] = 0
#normalize
norm = np.trapz(abs(psi_t)**2)
psi_t = psi_t/norm

#figure out the indeces at which we get to 1e-15s and 1e-16s to pause and save fig
time1 = 1e-15
time2 = 1e-16
ind1 = time1/h
ind1 = int(ind1)
ind2 = time2/h
ind2 = int(ind2)

figure()
for i in range(steps):
    #calculate v = B.psi(x,t)
    v = b1*psi_t[1:N] + b2*(psi_t[2:N+1] + psi_t[0:N-1])
    #solve for x using Gauss elim and backsub for tridiagnonal matrix
    psi_t[1:N] = banded(A, v, 1, 1)
    the_time = (i+1)*h
    #create animation
    clf()
    plot(x,psi_t.real)
    xlabel('$x$ meters')
    ylabel('$\psi(x,t)$')
    title('Wavefunction at time t = {0:.2e}'.format(the_time))
    draw()
    if i == ind1 or i == ind2:
        pause(10) #allow me to take screen shot
    else:
        pause(0.01)
    