#!/usr/bin/env python
# coding: utf-8

# In[11]:


import numpy as np
import astropy.units as u


# In[12]:


def Read(filename):
    
    #defining the function that takes the file as input
    
    file = open(filename, 'r')
    
    line1 = file.readline() #reads first line which is time
    label, value = line1.split()
    time = float(value) * u.Myr #changes to units of  10 Myear
    
    line2 = file.readline() #read second line which is number of particles
    label, value2 = line2.split()
    particles = float(value2)
   
    
    file.close()
    
    data = np.genfromtxt(filename, dtype= None, names=True, skip_header=3) #this stores the rest of the data file.
    
    #"dtype=None line is split using white spaces"
    #skip_header=3 skips the first three lines
    #names=true create arrays to store the data with the right labels
    #lables are "type, m, x, y, z, vx, vy, vz"

    
    return time, particles, data

   


# In[13]:


time, particles, data = Read('MW_000.txt') #this calls the function and tells it what to do 


# In[15]:


#print(data['vx'][1])


# In[5]:


#print(data['x'][1:])


# In[6]:


#print(particles)


# In[ ]:





# In[ ]:



