import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
sns.set()
# This program simulates Brownian motion in the presence of walls
# Note that the physical behaviour would be to stick to walls,
# which is the purpose of Q1a.
# Author: Nico Grisouard, University of Toronto
# Date: 14 November 2018

#################################################################
# This program simulates Brownian motion in the presence of walls
# Note that the physical behaviour would be to stick to walls,
# which is the purpose of Q1a.
# Author: Nico Grisouard, University of Toronto
# Date: 14 November 2018
#################################################################

def nextmove(x, y):
    """ randomly choose a direction
    1 = up, 2 = down, 3 = left, 4 = right"""
    direction =  np.random.randint(1,5)
    if direction == 1:  # move up
        y += 1
    elif direction == 2:  # move down
        y -= 1
    elif direction == 3:  # move right
        x += 1
    elif direction == 4:  # move left
        x -= 1
    else:
        print("error: direction isn't 1-4")
    if x == 100 :
        x -= 1
    elif y == 100 :
        y -= 1
    elif x == 0:
        x += 1
    elif y == 0:
        y += 1
    return x, y
    
Lp = 101  # size of domain
Nt = 5000  # number of time steps
# arrays to record the trajectory of the particle
x_position = np.empty(Nt)
y_position = np.empty(Nt)

centre_point = (Lp-1)//2  # middle point of domain
x_position[0] = centre_point
y_position[0] = centre_point
for i in range(Nt-1):
    x_position[i+1],y_position[i+1] = nextmove(x_position[i],y_position[i])
    
from pylab import clf, plot, xlim, ylim, show, pause
plt.figure(figsize = (10,10))  

for i in range(Nt):
    clf() # clear the plot
    plt.plot(x_position[0:i], y_position[0:i],c = 'brown')
    plt.ylim(0,100)
    plt.xlim(0,100)
    plt.draw()
    pause(0.01)
    #pause to allow a smooth animation
def nextmove(x, y):
    """ randomly choose a direction
    1 = up, 2 = down, 3 = left, 4 = right"""
    direction =  np.random.randint(1,5)
    if direction == 1:  # move up
        y += 1
    elif direction == 2:  # move down
        y -= 1
    elif direction == 3:  # move right
        x += 1
    elif direction == 4:  # move left
        x -= 1
    else:
        print("error: direction isn't 1-4")
    return x, y
x_particle_list = []
y_particle_list = []
Lp = 101  # size of domain
# arrays to record the trajectory of the particle
centre_point = (Lp-1)//2  # middle point of domain
x_position = np.array([centre_point])
y_position = np.array([centre_point])
j = 0
i = 0
while i < 100:   # number of particle count
    new_particle = False
    x_fill,y_fill = nextmove(x_position[j],y_position[j])
    if x_fill == Lp-1 or y_fill == Lp-1 or x_fill == 0 or y_fill == 0:
       # test for hitting a wall
        new_particle = True  
    for k in range(len(x_particle_list)):
        particle_position = np.array([x_fill,y_fill])
        end_x_array = x_particle_list[k]
        end_y_array = y_particle_list[k]
        rest_position = np.array([end_x_array[-1], end_y_array[-1]])
        if np.array_equal(particle_position, rest_position) == True:
            new_particle = True
    j += 1
    if new_particle == True:
       # print('generated new particle')
        # generate a new particle
        x_particle_list.append(x_position)
        y_particle_list.append(y_position)
        x_position = np.array([centre_point])
        y_position = np.array([centre_point])
        #print(len(x_position))
        i+= 1
        j = 0
    else:  
        x_position = np.append(x_position,x_fill)
        y_position = np.append(y_position,y_fill)
    # AND OFF YOU GO!
    # AND OFF YOU GO!
from pylab import clf, plot, xlim, ylim, show, pause
plt.figure(figsize = (10,10))  
for j in range(len(x_particle_list)):
    x_position = x_particle_list[j]
    y_position = y_particle_list[j]
    for i in range(len(x_position)):
        clf() # clear the plot
        if j > 0:
            clf() # clear the plot
            for k in range(j):
                x_rest = x_particle_list[k]
                y_rest = y_particle_list[k]
                plt.scatter(x_rest[-1],y_rest[-1],marker = "2" ,s = 80)
        plt.scatter(x_position[i], y_position[i*4],marker = "2" ,s = 80)
        plt.ylim(0,Lp-1)
        plt.xlim(0,Lp-1)
        plt.draw()
        pause(0.01)

#        
x_particle_list =[]
y_particle_list =[]
Lp = 151  # size of domain
# arrays to record the trajectory of the particle
centre_point = (Lp-1)//2  # middle point of domain
x_position = np.array([centre_point])
y_position = np.array([centre_point])
centre_full = False
new_particle = False


#
while centre_full == False:  # testing if centre point is full
    new_particle = False
    x_fill,y_fill = nextmove(x_position[j],y_position[j])
    if x_fill == Lp-1 or y_fill == Lp-1 or x_fill == 0 or y_fill == 0:
    # Here you test if particle hits a wall
        new_particle = True  
    for k in range(len(x_particle_list)):
        particle_position = np.array([x_fill,y_fill])
        end_x_array = x_particle_list[k]
        end_y_array = y_particle_list[k]
        rest_position = np.array([end_x_array[-1], end_y_array[-1]])
        if np.array_equal(particle_position, rest_position) == True:
            new_particle = True
        if np.array_equal(particle_position, np.array([centre_point,centre_point])) == True:
                centre_full = True
            # finishing the calculation if hit centre
    j += 1
    if new_particle == True:
       # print('generated new particle')
        # generate a new particle
        x_particle_list.append(x_position)
        y_particle_list.append(y_position)
        x_position = np.array([centre_point])
        y_position = np.array([centre_point])
        #print(len(x_position))
        i+= 1
        j = 0
    else:  
        x_position = np.append(x_position,x_fill)
        y_position = np.append(y_position,y_fill)
    # AND OFF YOU GO!
plt.figure(figsize = (10,10))  
for j in range(len(x_particle_list)):
    x_position = x_particle_list[j]
    y_position = y_particle_list[j]
    for i in range(len(x_position)):
        clf() # clear the plot
        if j > 0:
            clf() # clear the plot
            for k in range(j):
                x_rest = x_particle_list[k]
                y_rest = y_particle_list[k]
                plt.scatter(x_rest[-1],y_rest[-1],marker = "2" ,s = 80)
        plt.scatter(x_position[i], y_position[i*4],marker = "2" ,s = 80)
        plt.ylim(0,Lp-1)
        plt.xlim(0,Lp-1)
        plt.draw()
        pause(0.01)