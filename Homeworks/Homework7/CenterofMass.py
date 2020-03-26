#!/usr/bin/env python
# coding: utf-8

# In[12]:


# Homework 4
# Center of Mass Position and Velocity
# Jimena Stephenson


# ### Keep in mind this is just a template; you don't need to follow every step and feel free to change anything.
# ### We also strongly encourage you to try to develop your own method to solve the homework.

# In[13]:


# import modules
import numpy as np
import astropy.units as u
import astropy.table as tbl

from ReadFile import Read


# In[14]:


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
            
            xNew = self.x - XCOM2
            yNew = self.y - YCOM2
            zNew = self.z - ZCOM2
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
        
        COMPR = np.around(COMP, decimals = 2) * u.kpc
       
        return COMPR 


    def COM_V(self, COMX,COMY,COMZ):
        # Center of Mass velocity
        # input: X, Y, Z positions of the COM
        # returns 3D Vector of COM Velocities
        
        # the max distance from the center that we will use to determine the center of mass velocity                   
        RVMAX = 15.0 *u.kpc

        # determine the position of all particles relative to the center of mass position
        # write your own code below
        xV = self.x[:]*u.kpc - COMX
        yV = self.y[:]*u.kpc - COMY
        zV = self.z[:]*u.kpc - COMZ
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
            
        VXCOM, VYCOM, VZCOM = self.COMdefine(vxnew, vynew, vznew, mnew)

        # create a vector to store the COM velocity
        # set the correct units usint astropy
        # round all values
        # write your own code below
        COMV = [VXCOM, VYCOM, VZCOM]
        
        COMVR = np.round(COMV, decimals = 2) * u.km/u.s

        # return the COM vector                                                                                        
        return COMVR
    


# In[15]:


# Create a Center of mass object for the MW, M31 and M33
# below is an example of using the class for MW
MWCOM = CenterOfMass("MW_000.txt", 2)
M31COM = CenterOfMass("M31_000.txt", 2)
M33COM = CenterOfMass("M33_000.txt", 2)


# In[16]:


# MW:   store the position and velocity COM                                                                            
MW_COM_P = MWCOM.COM_P(0.1)
MW_COM_V = MWCOM.COM_V(MW_COM_P[0],MW_COM_P[1],MW_COM_P[2])

# M31:   store the position and velocity COM                                                                           
M31_COM_P = M31COM.COM_P(0.1)
M31_COM_V = M31COM.COM_V(M31_COM_P[0],M31_COM_P[1],M31_COM_P[2])


# M33:   store the position and velocity COM                                                                           
M33_COM_P = M33COM.COM_P(0.1)
M33_COM_V = M33COM.COM_V(M33_COM_P[0],M33_COM_P[1],M33_COM_P[2])


# In[17]:


# Making the Table of Results

# we first gather all results and put them into a numpy array
# *tuple() will extract the tuple, kind of like C-language
tab_results = ['MW COM', *tuple(MW_COM_P.value), *tuple(MW_COM_V.value), 
               'M31 COM', *tuple(M31_COM_P.value), *tuple(M31_COM_V.value), 
               'M33 COM', *tuple(M33_COM_P.value), *tuple(M33_COM_V.value)]
tab_results = np.reshape(tab_results, (3, 7)) # convert 1D array to 2D

# make and display the table with astropy
t = tbl.Table(tab_results, 
              names = ['Galaxy', 'X-component (kpc)', 'Y-Component (kpc)', 
                       'Z-Component (kpc)','VX-component (km/s)', 
                       'VY-Component (km/s)', 'VZ-Component (km/s)'])
t.show_in_notebook()


# In[18]:


# Q2 
# Determine the separation between the MW and M31                                                                      
MW_M31 = np.sqrt((M31_COM_P[0]-MW_COM_P[0])**2 + (M31_COM_P[1]-MW_COM_P[1])**2 + (M31_COM_P[2]-MW_COM_P[2])**2)
print("Separation between the MW and M31 =", np.round(MW_M31))

# Determine the relative velocity between the MW and M31                                                                      
vMW_M31 = np.sqrt((M31_COM_V[0]-MW_COM_V[0])**2 + (M31_COM_V[1]-MW_COM_V[1])**2 + (M31_COM_V[2]-MW_COM_V[2])**2)
print("Relative Velocity between the MW and M31 =", np.round(vMW_M31))


# In[19]:


# Q3
# Determine the relative position between M33 and M31                                                                  
M33_M31 = np.sqrt((M33_COM_P[0]-M31_COM_P[0])**2 + (M33_COM_P[1]-M31_COM_P[1])**2 + (M33_COM_P[2]-M31_COM_P[2])**2)
print("Relative Position between M33 and M31 = ", np.round(M33_M31))


# Determine the relative velocity between M33 and M31                                                                  
vM33_M31 = np.sqrt((M33_COM_V[0]-M31_COM_V[0])**2 + (M33_COM_V[1]-M31_COM_V[1])**2 + (M33_COM_V[2]-M31_COM_V[2])**2)
print("Relative Velocity between M33 and M31 = ", np.round(vM33_M31))


# Q4: The iterative procedue is necessary because as galaxies interact their stars get thrown to large radii with different speeds. Tidal tails are not spherically symmetric; this can cause the COM position calculation to be really off if you are using all the particles at large radii.

# In[ ]:



