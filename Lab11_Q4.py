# -*- coding: utf-8 -*-
"""
Created on Wed Nov 28 11:56:54 2018

@author: brayd
"""
#q4
#Brayden Kell
#   Adapted Prof Grisouard's code salesman_start.py
#       Nov 28 2018
import numpy as np
import matplotlib.pyplot as plt
from random import random, randrange, seed
#from pdb import set_trace
# %% Main program starts here ------------------------------------------------|
N = 50
Tmax = 10.0
Tmin = 1e-3
tau = 1e3 #varied this 

# Choose N city locations and calculate the initial distance
dimer_dic = dict()
#set up 50x50 grid space
X = int(N)
Y = int(N)
'''
# Set up the graphics
plt.figure(1)
plt.plot(r[:, 0], r[:, 1], 'o-',  markersize=6, linewidth=2, color='k')
plt.draw()
plt.pause(0.001)
'''
# Main loop
t = 0
T = Tmax
Energy = 1
totEnergy_atTime = []
seed(91) #Vary this, get a few outputs
d = 0 #counter for indexing dictionary key
while T > Tmin:
    # Cooling
    t += 1
    T = Tmax*np.exp(-t/tau)
    '''
    # Update the visualization every 100 moves
    if t % 100 == 0:
        plt.figure(1)
        plt.clf()
        plt.plot(r[:, 0], r[:, 1], 'o-',  markersize=6, linewidth=2, color='k')
        plt.draw()
        plt.pause(0.001)
    '''
    #Choose two adjacent sites
    i1, j1 = randrange(1, X-1), randrange(1, Y-1) #first site in [1,49]x[1,49] to prevent going outside [50,50]x[50,50]
    i2, j2 = randrange(i1-1, i1+1), randrange(j1-1, j1+1) #choose adjacent grid coordinate
    while (i2 == i1 and j2 == j1) or not(i2 == i1 or j2 == j1 ): #make sure it's distinct and not diagonal
        i2, j2 = randrange(i1-1, i1+1), randrange(j1-1, j1+1)
    site1 = [i1,j1]
    site2 = [i2,j2]
        
    freeSpace = True #flag indicating space is free to add dimer
    toDelete = False #flag indicating a space is occupied and sample exceeded boltzman prop. so delete it
    for k in dimer_dic.keys():
        theDimer = dimer_dic[k]
        if [site1,site2] == theDimer or [site2,site1] == theDimer: #occupied by a single dimer
            freeSpace = False
            if random() < np.exp(-Energy/T): #if exceed probability, delete dimer
                toDelete = True
                #set_trace()
                ind = k
        elif theDimer[0] == site1 or theDimer[0] == site2 or theDimer[1] == site1 or theDimer[1] == site2:
            freeSpace = False
    if freeSpace: #if unoccupied add dimer
        d += 1 #add one to counter
        dimer_dic[d] = [site1,site2]
    if toDelete:#if occupied and passed boltzmann
        del dimer_dic[ind]
        
    totEnergy_atTime.append(len(dimer_dic)) #update energy for current time point

plt.figure()
for k in dimer_dic.keys():
    theDimer = dimer_dic[k]
    plt.plot([theDimer[0][0],theDimer[1][0]], [theDimer[0][1],theDimer[1][1]],'o-',\
                 markersize=2, linewidth=2)#plot final configuration 
    
x = np.arange(len(totEnergy_atTime))
plt.figure()
plt.plot(x, totEnergy_atTime)
plt.title('Dimers vs. Time')
plt.xlabel('Time (1/iteration)')
plt.ylabel('# of dimers')
plt.axhline(y=50*50/2, color='r', linestyle='-')

print('Number of dimers in lattice is ',len(dimer_dic))
