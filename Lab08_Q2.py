import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import seaborn as sns
sns.set()
# Defining Constants
h = 10e-6  # time step
L = 1   # Length of wire in metres
v = 100 # velocity in m/s
d = 0.1 # distance of hammer
C = 1   # 1 m/s
sigma = 0.3  # sigma in metres
N = 100 # grid spacings
a = L/N

dt= 0.003


dt= dt= 0.01
iterations = int(dt/h)
displacement = np.zeros([iterations,N+1],float)
velocity = np.zeros([iterations,N+1],float)
x = np.linspace(0,L,N+1)
velocity[0,:] = C * x*(L-x)/L**2*np.exp( -1*(x-d)**2/(2*sigma**2))

for i in range(iterations-1): # time iteration
    for j in range(1,len(x)-1): # space iterations
        displacement[i+1,j] = displacement[i,j] +  h*velocity[i,j]
        velocity[i+1,j] = velocity[i,j] + \
        h*(v**2)/(a**2)*(displacement[i,j+1] + displacement[i,j-1] - 2*displacement[i,j])  


from pylab import clf, plot, xlim, ylim, show, pause
'''
for i in range(iterations):
	clf() # clear the plot
	plt.plot(displacement[i*10,:])
	ylim([-1, 1]) # set the x boundaries constant
	#xlim([0, 500 ]) # and the y boundaries
	plt.draw()
	pause(0.01)
	#pause to allow a smooth animation
'''
from pylab import clf, plot, xlim, ylim, show, pause
for i in range((iterations//10)-1):
	clf() # clear the plot
	plt.plot(x,velocity[i*10,:])#,s = 3)
	ylim([-.2, .2]) # set the x boundaries constant
	#xlim([0, 500 ]) # and the y boundaries
	plt.draw()
	plt.title('Velocity Time Evolution', fontsize = 12)
	plt.xlabel('Position (m)', fontsize = 10)
	plt.ylabel('Velocity (m/s)',  fontsize = 10)
	pause(0.01)
	#pause to allow a smooth animation

from pylab import clf, plot, xlim, ylim, show, pause
for i in range((iterations//10)-1):
	clf() # clear the plot
	plt.plot(x,displacement[i*10,:])#,s = 3)
	ylim([-.0015, .0015]) # set the x boundaries constant
	#xlim([0, 500 ]) # and the y boundaries
	plt.title('Displcement Time Evolution',fontsize = 12)
	plt.xlabel('Position (m)', fontsize = 10)
	plt.ylabel('displacement (m)', fontsize = 10)
	plt.draw()
	pause(0.01)
	#pause to allow a smooth animation
