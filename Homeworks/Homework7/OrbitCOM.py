#!/usr/bin/env python
# coding: utf-8

# In[10]:


# Homework 6 Template
# Jimena Stephenson 


# In[1]:


# import modules
import numpy as np
import astropy.units as u
from astropy.constants import G

# import plotting modules
import matplotlib.pyplot as plt
import matplotlib
get_ipython().run_line_magic('matplotlib', 'inline')

# my modules
from ReadFile import Read
from CenterOfMass2 import CenterOfMass


# In[4]:


def OrbitCom(galaxy, start, end, n=5):
    """function that loops over all the desired snapshots to compute the COM pos and vel as a function of time.
    inputs:
    galaxy name ie: "MW", start: first snapshot to be read, end: last snapshot to be read, n: intervals to return the COM
    We will compute up to snapshot 800 and we will outpu values in intervals of n=5      
    returns: 
    computes the timen and COM position and velocity vectors of a given a galaxy in each snapshot
    saves that output into a file.
    """
    
    # compose the filename for output
    fileout = "Orbit_"+galaxy+".txt"
    
    #  set tolerance and VolDec for calculating COM_P in CenterOfMass
    # for M33 that is stripped more, use different values for VolDec
    
    if (galaxy == 'M33'):
        delta = 0.1
        VolDec=4
        
    else:
        delta = 0.1
        VolDec = 2

    
    # generate the snapshot id sequence 
    # it is always a good idea to also check if the input is eligible (not required)
    snap_ids = np.arange(start, end+1, step=n)
    
    
    # initialize the array for orbital info: t, x, y, z, vx, vy, vz of COM
    
    orbit = np.zeros((snap_ids.size, 7))
    
    # a for loop 
    for  i, snap_id in enumerate(snap_ids):
        
        # compose the data filename (be careful about the folder)
         # inputs to recontruct the filenumber to the value "000"
        ilbl = '000' + str(snap_id)
        # remove all but the last 3 digits
        ilbl = ilbl[-3:]
        filename = "%s_"%(galaxy) + ilbl + '.txt'
        
        # Initialize an instance of CenterOfMass class, using disk particles
        
        COM = CenterOfMass(filename, 2)

        # Store the COM pos and vel. Remember that now COM_P required VolDec
        COM_P = COM.COM_P(delta, VolDec)
        COM_V = COM.COM_V(COM_P[0],COM_P[1],COM_P[2])

    
        # store the time, pos, vel in ith element of the orbit array,  without units (.value) 
        # note that you can store 
        # a[i] = var1, *tuple(array1)
        
        orbit[i] = COM.time.value/1000, *tuple(COM_P.value), *tuple(COM_V.value)
        
        # print snap_id to see the progress
        print(snap_id)
        
    # write the data to a file
    # we do this because we don't want to have to repeat this process 
    # this code should only have to be called once per galaxy.
    np.savetxt(fileout, orbit, fmt = "%11.3f"*7, comments='#',
               header="{:>10s}{:>11s}{:>11s}{:>11s}{:>11s}{:>11s}{:>11s}"\
                      .format('t', 'x', 'y', 'z', 'vx', 'vy', 'vz'))


# In[12]:


OrbitCom('MW', 0, 801, n=5)
OrbitCom('M31', 0, 801, n=5)
OrbitCom('M33', 0, 801, n=5)


# In[20]:


# Read in the data files for the orbits of each galaxy that you just created
# headers:  t, x, y, z, vx, vy, vz
# using np.genfromtxt
OrbitMW_data = np.genfromtxt("Orbit_MW.txt")
OrbitM31_data = np.genfromtxt("Orbit_M31.txt")
OrbitM33_data = np.genfromtxt("Orbit_M33.txt")


# In[17]:


# function to compute the magnitude of the difference between two vectors 
# You can use this function to return both the relative position and relative velocity for two 
# galaxies over the entire orbit  


def vect_diff(vc1, vc2):
    
    diff = vc1-vc2
    
    return np.sqrt(diff[:,0]**2+diff[:,1]**2+diff[:,2]**2)


# In[18]:


# Determine the magnitude of the relative position and velocities 

# of MW and M31
MWM31 = vect_diff(OrbitMW_data[:, 1:4], OrbitM31_data[:, 1:4])
VMWM31 = vect_diff(OrbitMW_data[:, 4:8], OrbitM31_data[:, 4:8])

# of M33 and M31
M33M31 = vect_diff(OrbitM33_data[:, 1:4], OrbitM31_data[:, 1:4])
VM33M31 = vect_diff(OrbitM33_data[:, 4:8], OrbitM31_data[:, 4:8])


# In[22]:


# Plot the Orbit of the galaxies 
#################################

### MW - M31 SEPARATION

fig = plt.figure(figsize=(10,10))
ax = plt.subplot(111)
plt.plot(OrbitMW_data[:,0], MWM31, color='turquoise', linewidth=4, label='MW-M31 Sepparation')

# Add axis labels
plt.xlabel('Time (Gyr)', fontsize=22)
plt.ylabel('Separation (kpc)', fontsize=22)

#adjust tick label font size
label_size = 22
matplotlib.rcParams['xtick.labelsize'] = label_size 
matplotlib.rcParams['ytick.labelsize'] = label_size

# add a legend with some customizations.
legend = ax.legend(loc='upper right',fontsize='x-large')

## M33-M31 SEPARATION
fig = plt.figure(figsize=(10,10))
ax = plt.subplot(111)
plt.plot(OrbitM31_data[:,0], M33M31, color='teal', linewidth=4, label='M33-M31 Sepparation')

# Add axis labels
plt.xlabel('Time (Gyr)', fontsize=22)
plt.ylabel('Separation (kpc)', fontsize=22)


#adjust tick label font size
label_size = 22
matplotlib.rcParams['xtick.labelsize'] = label_size 
matplotlib.rcParams['ytick.labelsize'] = label_size

# add a legend with some customizations.
legend = ax.legend(loc='upper right',fontsize='x-large')



# In[24]:


# Plot the orbital velocities of the galaxies 
#################################

### MW - M31 ORBITAL VEL

fig = plt.figure(figsize=(10,10))
ax = plt.subplot(111)
plt.plot(OrbitMW_data[:,0], VMWM31, color='salmon', linewidth=4, label='MW-M31 Orbital Vel')

# Add axis labels
plt.xlabel('Time (Gyr)', fontsize=22)
plt.ylabel('Relative Velocity (km/s)', fontsize=22)


#adjust tick label font size
label_size = 22
matplotlib.rcParams['xtick.labelsize'] = label_size 
matplotlib.rcParams['ytick.labelsize'] = label_size

# add a legend with some customizations.
legend = ax.legend(loc='upper right',fontsize='x-large')

### M33 - M31 ORBITAL VELOCITY

fig = plt.figure(figsize=(10,10))
ax = plt.subplot(111)
plt.plot(OrbitM31_data[:,0], VM33M31, color='plum', linewidth=4, label='M33-M31 Orbital Vel')

# Add axis labels
plt.xlabel('Time (Gyr)', fontsize=22)
plt.ylabel('Relative Velocity (km/s)', fontsize=22)


#adjust tick label font size
label_size = 22
matplotlib.rcParams['xtick.labelsize'] = label_size 
matplotlib.rcParams['ytick.labelsize'] = label_size

# add a legend with some customizations.
legend = ax.legend(loc='upper right',fontsize='x-large')


# QUESTIONS

# Q1 It seems that MW and M31 will have to close encounters before merging

# Q2 The relative velocity dies down as the merge

# Q3 MW and M31 merge in about 6.2 Gyr. M33 orbits closer when the other two merge

# In[ ]:




