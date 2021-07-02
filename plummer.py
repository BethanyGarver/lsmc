import numpy as np

# This program returns a plummer model of N stars
# with the given mass and spatial scale.  The velocity
# scale is set by the mass and spatial scale.
# 
# If rscale is not set then a scale is computed for
# an average "fluffy" galaxy us the relation:
#
#   rscale = 0.9 * (Mass/1e9)^(1/3)
#
# INPUT
#
#   N       Number of particles
#   MASS    Total mass of system in solar masses
#   RSCALE  Spatial scale of galaxy in kpc
#   /SOFT   Put in a column for softening parameters (set to 0.)
#
# OUTPUT
#
#   ARR     Array of particles and their properties, [N,7 or 8]
#           soft = 0,   mass, x, y, z, vx, vy, vz
#           soft = 1,   mass, 0, x, y, z, vx, vy, vz

def plummer(n=100000,mass=1.,rscale=None,soft=0,xyz=[0,0,0],v=[0,0,0]):

    if rscale is None:
        rscale = 0.9*(mass/1e9)**(1./3.)
    
    xarr = np.tile(xyz,(n,1))
    varr = np.tile(v,(n,1))

    # getting the data from "BI"
    biarr = np.loadtxt('/home/bgarver/Documents/galsat/data/plummer.data',skiprows=1)

    # put in a softening parameter column
    if soft==1:

        arr = np.zeros((n,8))
        arr[:,0] = mass/float(n)
        arr[:,2:5] = biarr[0:n,1:4]  #x,y,z
        arr[:,5:8] = biarr[0:n,4:7]  #vx,vy,vz

        # scaling the spatial part
        rold = 1.5
        arr[:,2:5] = arr[:,2:5]*rscale/rold+xarr

        # scaling the velocity part
        G = 4.4967e-6              # G in kpc/Gyr^2
        vscale = np.sqrt(G*mass/rscale)
        arr[:,5:8] = arr[:,5:8]*vscale+varr

    # no softening parameter column
    else:

        arr = np.zeros((n,7))
        arr[0,:] = mass/float(n)
        arr[:,1:4] = biarr[0:n,1:4]   #x,y,z
        arr[:,4:7] = biarr[0:n,4:7]   #vx,vy,vz

        # scaling the spatial part
        rold = 1.5
        arr[:,1:4] = arr[:,1:4]*rscale/rold+xarr

        # scaling the velocity part
        G = 4.4967e-6              # G in kpc/Gyr^2
        vscale = np.sqrt(G*mass/rscale)
        arr[:,4:7] = arr[:,4:7]*vscale+varr
    
    return arr