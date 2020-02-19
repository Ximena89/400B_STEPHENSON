#!/usr/bin/env python
# coding: utf-8

# In[5]:


# Homework 4
# Center of Mass Position and Velocity
# Jimena Stephenson


# ### Keep in mind this is just a template; you don't need to follow every step and feel free to change anything.
# ### We also strongly encourage you to try to develop your own method to solve the homework.

# In[6]:


# import modules
import numpy as np
import astropy.units as u
import astropy.table as tbl

from ReadFile import Read


# In[110]:


class CenterOfMass:
# Class to define COM position and velocity properties of a given galaxy 
# and simulation snapshot
    
    
    def __init__(self, filename, ptype):
    # Initialize the instance of this Class with the following properties:
    
        # read data in the given file using Read
        self.time, self.total, self.data = Read(filename)                                                                                             

        #create an array to store indexes of particles of desired Ptype                                
        self.index = np.where(self.data['type'] == ptype)

        # store the mass, positions, velocities of only the particles of the given type
        # the following only gives the example of storing the mass
        self.m = self.data['m'][self.index]
        # write your own code to complete this for positions and velocities
       
        #Positions
        self.x = self.data['x'][self.index]
        self.y = self.data['y'][self.index]
        self.z = self.data['z'][self.index]
        
        #Velocities
        self.vx = self.data['vx'][self.index]
        self.vy = self.data['vy'][self.index]
        self.vz = self.data['vz'][self.index]

    def COMdefine(self,a,b,c,m):
    # Function to compute the center of mass position or velocity generically
    # input: array (a,b,c) of positions or velocities and the mass
    # returns: 3 floats  (the center of mass coordinates)

        # write your own code to compute the generic COM using Eq. 1 in the homework instructions
        
        # xcomponent Center of mass
        na = np.sum(a*m)
        da = np.sum(m)
        Acom = na/da
        
        # ycomponent Center of mass
        nb = np.sum(b*m)
        db = np.sum(m)
        Bcom = nb/db 
        
        # zcomponent Center of mass
        nc = np.sum(c*m)
        dc = np.sum(m)
        Ccom = nc/dc
        
        return Acom, Bcom, Ccom
    
 

   
    def COM_P(self, delta):
    # Function to specifically return the center of mass position and velocity                                         
    # input:                                                                                                           
    #        particle type (1,2,3)                                                                                     
    #        delta (tolerance)                                                                                         
    # returns: One vector, with rows indicating:                                                                                                                                                                            
    #       3D coordinates of the center of mass position (kpc)                                                             

        # Center of Mass Position                                                                                      
        ###########################                                                                                    

        # Try a first guess at the COM position by calling COMdefine                                                   
        XCOM, YCOM, ZCOM = self.COMdefine(self.x, self.y, self.z, self.m)
        # compute the magnitude of the COM position vector.
        # write your own code below
        RCOM = np.sqrt(XCOM**2 + YCOM**2 + ZCOM**2)


        # iterative process to determine the center of mass                                                            

        # change reference frame to COM frame                                                                          
        # compute the difference between particle coordinates                                                          
        # and the first guess at COM position
        # write your own code below
        xNew = self.x - XCOM
        yNew = self.y - YCOM
        zNew = self.z - ZCOM
        RNEW = np.sqrt(xNew**2 + yNew**2 + zNew**2)

        # find the max 3D distance of all particles from the guessed COM                                               
        # will re-start at half that radius (reduced radius)                                                           
        RMAX = max(RNEW)/2.0
        
        # pick an initial value for the change in COM position                                                      
        # between the first guess above and the new one computed from half that volume
        # it should be larger than the input tolerance (delta) initially
        CHANGE = 1000.0 #kpc

        # start iterative process to determine center of mass position                                                 
        # delta is the tolerance for the difference in the old COM and the new one.    
        
        while (CHANGE > delta):
            # select all particles within the reduced radius (starting from original x,y,z, m)
            # write your own code below (hints, use np.where)
            index2 = np.where(RNEW <= RMAX)
            x2 = self.x[index2]
            y2 = self.y[index2]
            z2 = self.z[index2]
            m2 = self.m[index2]

            # Refined COM position:                                                                                    
            # compute the center of mass position using                                                                
            # the particles in the reduced radius
            # write your own code below

        
            XCOM2, YCOM2, ZCOM2 = self.COMdefine(x2, y2, z2, m2)
            
            
            # compute the new 3D COM position
            # write your own code below
            RCOM2 = np.sqrt(XCOM2**2 + YCOM2**2 + ZCOM2**2)

            # determine the difference between the previous center of mass position                                    
            # and the new one.                                                                                         
            CHANGE = np.abs(RCOM - RCOM2)
            # uncomment the following line if you wnat to check this                                                                                               
            #print ("CHANGE = ", CHANGE)                                                                                     

            # Before loop continues, reset : RMAX, particle separations and COM                                        

            # reduce the volume by a factor of 2 again                                                                 
            RMAX = RMAX/2.0
            # check this. 
            #print ("maxR", RMAX)

            # Change the frame of reference to the newly computed COM.                                                 
            # subtract the new COM
            # write your own code below
            
            xNew = x2 - XCOM2
            yNew = y2 - YCOM2
            zNew = z2 - ZCOM2
            RNEW = np.sqrt(xNew**2 + yNew**2 + zNew**2)


            # set the center of mass positions to the refined values                                                   
            XCOM = XCOM2
            YCOM = YCOM2
            ZCOM = ZCOM2
            RCOM = RCOM2

            # create a vector to store the COM position                                                                                                                                                       
            COMP = [XCOM, YCOM, ZCOM]

        # set the correct units using astropy and round all values
        # and then return the COM positon vector
        # write your own code below
        
        COMPR = np.around(COMP, decimals = 2) 
       
       
            
        return COMPR 


    def COM_V(self, COMX,COMY,COMZ):
        # Center of Mass velocity
        # input: X, Y, Z positions of the COM
        # returns 3D Vector of COM Velocities
        
        # the max distance from the center that we will use to determine the center of mass velocity                   
        RVMAX = 15.0 #*u.kpc

        # determine the position of all particles relative to the center of mass position
        # write your own code below
        xV = self.x - COMX
        yV = self.y - COMY
        zV = self.z - COMZ
        RV = np.sqrt(xV**2 + yV**2 +zV**2)
        
        # determine the index for those particles within the max radius
        # write your own code below
        indexV = np.where(RV < RVMAX)

        # determine the velocity and mass of those particles within the max radius
        # write your own code below
        vxnew = self.vx[indexV]
        vynew = self.vy[indexV]
        vznew = self.vz[indexV]
        mnew =  self.m[indexV] 
        
        # compute the center of mass velocity using those particles
        # write your own code below
        
         # xcomponent Center of mass
        nav = np.sum(vxnew*mnew)
        dav = np.sum(mnew)
        Acomv = nav/dav
        
            # ycomponent Center of mass
        nbv = np.sum(vynew*mnew)
        dbv = np.sum(mnew)
        Bcomv = nbv/dbv 
        
            # zcomponent Center of mass
        ncv = np.sum(vznew*mnew)
        dcv = np.sum(mnew)
        Ccomv = ncv/dcv
            
        VXCOM, VYCOM, VZCOM = Acomv, Bcomv, Ccomv

        # create a vector to store the COM velocity
        # set the correct units usint astropy
        # round all values
        # write your own code below
        COMV = [VXCOM, VYCOM, VZCOM]
        
        COMVR = np.round(COMV, decimals = 2)

        # return the COM vector                                                                                        
        return COMVR
    


# In[111]:


# Create a Center of mass object for the MW, M31 and M33
# below is an example of using the class for MW
MWCOM = CenterOfMass("MW_000.txt", 2)
M31COM = CenterOfMass("M31_000.txt", 2)
M33COM = CenterOfMass("M33_000.txt", 2)


# In[ ]:





# In[112]:


# below gives you an example of calling the class's functions
# MW:   store the position and velocity COM
MW_COMP = MWCOM.COM_P(0.1)
MW_COMV = MWCOM.COM_V(MW_COMP[0], MW_COMP[1], MW_COMP[2])
print(MW_COMP)
print(MW_COMV)


# Question #1
# 
# The center of mass of the MW is close to the origin

# In[113]:


# now write your own code to answer questions


# Center of mass and velocity for M31

# In[114]:


M31_COMP = M31COM.COM_P(0.1)
M31_COMV = M31COM.COM_V(M31_COMP[0], M31_COMP[1], M31_COMP[2])
print(M31_COMP)
print(M31_COMV)


# Center of mass for M33

# In[115]:


M33_COMP = M33COM.COM_P(0.1)
M33_COMV = MWCOM.COM_V(MW_COMP[0], MW_COMP[1], MW_COMP[2])
print(M33_COMP)
print(M33_COMV)


# Question #2
# 
# Magnitude of the current separation between MW and M31

# In[123]:


#Calculating the separation 
magMW = np.sqrt(MW_COMP[0]**2 + MW_COMP[1]**2 + MW_COMP[2]**2)
magM31 = np.sqrt(M31_COMP[0]**2 + M31_COMP[1]**2 + M31_COMP[2]**2)

#Calculating the velocity
magMWv = np.sqrt(MW_COMV[0]**2 + MW_COMV[1]**2 + MW_COMV[2]**2)
magM31v = np.sqrt(M31_COMV[0]**2 + M31_COMV[1]**2 + M31_COMV[2]**2)

separation = magMW - magM31 
velocity = magMWv - magM31v

print(separation)
print(velocity)


# Question #3
# 
# Magnitude of separation between M31 and M33

# In[132]:


#Calculating the separation 
magM31 = np.sqrt(M31_COMP[0]**2 + M31_COMP[1]**2 + M31_COMP[2]**2)
magM33 = np.sqrt(M33_COMP[0]**2 + M33_COMP[1]**2 + M33_COMP[2]**2)

#Calculating the velocity
magM31v = np.sqrt(M31_COMV[0]**2 + M31_COMV[1]**2 + M31_COMV[2]**2)
magM33v = np.sqrt(M33_COMV[0]**2 + M33_COMV[1]**2 + M33_COMV[2]**2)

separation = magM33 - magM31 
velocity = magM33v - magM31v

print(separation)
print(velocity)


# I had a little truble with the units. Jupiter didn't like at any point where I tried to use astropy units. It would not let me do math with units.

# Question #4
# 
# It is important beacuse as the move and interact gravitationally with each other, the center of mass keeps changing

# In[ ]:




