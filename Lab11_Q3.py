# Dimer covering problem (Newman's 10.11) solution
# Author: Nico Grisouard, University of Toronto
# Adapted from Newman's salesman.py
# Date: 22 November 2018

#q3
#Brayden Kell
#   Adapted Prof Grisouard's code
#       Nov 28 2018

import numpy as np
import matplotlib.pyplot as plt
from random import random, randrange, seed


def mag(x):
    """ Function to calculate the magnitude of a vector """
    return np.sqrt(x[0]**2 + x[1]**2)


def distance():
    """ Function to calculate the total length of the tour """
    s = 0.0
    for i in range(N):
        s += mag(r[i+1] - r[i])
    return s


# %% Main program starts here ------------------------------------------------|
N = 25
R = 0.02
Tmax = 10.0
Tmin = 1e-3
tau = 1e4 #varied this 


# Choose N city locations and calculate the initial distance
r = np.empty([N+1, 2], float)
seed(1738) #for reproducibility 
for i in range(N):
    r[i, 0] = random()
    r[i, 1] = random()
r[N] = r[0]
D = distance()

# Set up the graphics
plt.figure(1)
plt.plot(r[:, 0], r[:, 1], 'o-',  markersize=6, linewidth=2, color='k')
plt.draw()
plt.pause(0.001)

# Main loop
t = 0
T = Tmax
seed(19) #Vary this, get a few outputs
while T > Tmin:

    # Cooling
    t += 1
    T = Tmax*np.exp(-t/tau)

    # Update the visualization every 100 moves
    if t % 100 == 0:
        plt.figure(1)
        plt.clf()
        plt.plot(r[:, 0], r[:, 1], 'o-',  markersize=6, linewidth=2, color='k')
        plt.draw()
        plt.pause(0.001)

    # Choose two cities to swap and make sure they are distinct
    i, j = randrange(1, N), randrange(1, N)
    while i == j:
        i, j = randrange(1, N), randrange(1, N)

    # Swap them and calculate the change in distance
    oldD = D
    r[i, 0], r[j, 0] = r[j, 0], r[i, 0]
    r[i, 1], r[j, 1] = r[j, 1], r[i, 1]
    D = distance()
    deltaD = D - oldD

    # If the move is rejected, swap them back again
    if random() > np.exp(-deltaD/T):
        r[i, 0], r[j, 0] = r[j, 0], r[i, 0]
        r[i, 1], r[j, 1] = r[j, 1], r[i, 1]
        D = oldD
   
print('Final distance is ', D) 
    
