# -*- coding: utf-8 -*-
"""
Created on Wed Oct 31 11:21:21 2018

@author: brayd
"""
'''
Using Gauss-Seidel method with overrelaxation to calculate steady
state heat distribution for a particular 2D geometric object with some particular
boundary conditions (see figure 1 in lab manual). We produce an animation 
for an illustration of the iterative process towards convergence within some
pre-defined error tolerance. 
'''
from numpy import zeros,max,ones, copy
from pylab import imshow, clf, pause, draw, transpose, gray
import matplotlib.pyplot as plt

# Constants
a = 0.1 #cm - grid spacing
M = int(20/a) # number of grid squares on horizontal side
N = int(8/a) # number of grid squares on vertical side 
x_gap = int(5/a) #length of segment AB and EF
y_gap = int(3/a) #length of BC and DE
target = 1e-6   # Target accuracy
w = 0.9 #overrelaxation parameter 

# Create arrays to hold potential values
phi = zeros([M+1,N+1],float) #hold temperature values at each grid coordinate
A = ones([M+1,N+1],int) #we will use this to store 0s in regions we dont
#want to update with Gauss-Seidel method and 1's where we do

#set boundary conditions
#notice we also set A equal to zero in the grids on the boundary and empty region
for i in range(M+1):
    for j in range(N+1):
        if i<=x_gap and j == 0: #AB
            phi[i,j] = i*a
            A[i,j] = 0
        elif i==x_gap and j<=y_gap: #BC
            phi[i,j] = 5 + (2/3)*j*a
            A[i,j] = 0
        elif i>=x_gap and i<=(M - x_gap) and j == y_gap: #CD
            phi[i,j] = 7
            A[i,j] = 0
        elif i==(M - x_gap) and j<=y_gap: #DE
            phi[i,j] = 5 + (2/3)*j*a
            A[i,j] = 0
        elif i>=(M - x_gap) and j == 0: #EF
            phi[i,j] = 20 - i*a
            A[i,j] = 0
        elif i == M: #FG
            phi[i,j] = (10/8)*j*a
            A[i,j] = 0
        elif j == N: #GH
            phi[i,j] = 10
            A[i,j] = 0
        elif i == 0: #HA
            phi[i,j] = (10/8)*j*a
            A[i,j] = 0
        if i>=x_gap and i<=(M - x_gap) and j < y_gap: #empty region
            A[i,j] = 0
# Main loop
delta = 1.0 #an initial delta to enter while loop
ind = 0 #counting number of iterations
while delta>target: #for part c
#while ind<100: #for part b
    phi0 = copy(phi) #copy the phi from the previous iterations to calculate
    #difference at end of iteration
    ind = ind+1
    # Calculate new values of the potential
    for i in range(M+1):
        for j in range(N+1):
            if A[i,j] == 1: 
                phi[i,j] = ((1+w)/4)*(phi[i+1,j] + phi[i-1,j] + phi[i, j+1] +\
                   phi[i, j-1]) - w*phi[i,j]
    # Calculate maximum difference from old values
    delta = max(abs(phi-phi0))

    #create animation
    clf()
    imshow(transpose(phi))
    gray()
    #flip y axis so the shape looks as in lab manual (doesnt change the physics)
    plt.gca().invert_yaxis() 
    draw()
    pause(0.01)

print('The number of iterations is ', ind) 
print('The temperature at x=2.5cm, y=1cm in celsuis is ', phi[25,10]) 