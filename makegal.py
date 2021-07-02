from plummer import plummer
import numpy as np

def makegal(input0,nbody,noprint=None,stp=None,dhalo=None):


    #+
    # This program runs an N-test-body simulation of a
    # galaxy with a plummer distribution of stars that
    # orbits the Milky Way.
    #
    # INPUTS:
    #  input    Input parameters [mass,soft,x,y,z,vx,vy,vz] of the galaxy CM
    #  =nbody   Number of bodies to use in N-test-body simulation
    #  /noprint Don't print out the parameters
    #  /stp     Stop at the end of the program
    #
    # OUTPUTS:
    #  arr      The array of parameters for the galaxy.
    #
    # CALLS:
    #  plummer   To get the distribution of the N-bodies
    #
    # By David Nidever   March 2006
    #-

    # Not enough inputs
    if len(input0)==0:
        print('Syntax - makegal,input,arr,nbody=nbody,noprint=noprint,stp=stp')
        return

    # INITIAL CONDITIONS
    #input = [mass,soft,xstart,ystart,zstart,vxstart,vystart,vzstart]
    mass = input0[0]
    soft = input0[1]
    x0 = input0[2]
    y0 = input0[3]
    z0 = input0[4]
    vx0 = input0[5]
    vy0 = input0[6]
    vz0 = input0[7]


    # GETTING THE PLUMMER DISTRIBUTION OF STARS
    #mass = 1e7
    #rscale = 0.9*(mass/1d9)^(1./3.)

    print('USING {} STARS'.format(nbody))
    print('Getting Plummer distribution of stars')
    arr = plummer(n=nbody,mass=mass,rscale=soft,soft=1)


    # VELOCITIES FROM PLUMMER ARE IN KPC/GYR
    # CONVERT TO KM/S FOR GALSAT INPUT

    # CONVERTING VELOCITIES (kpc/Gyr -> km/s)
    # 1 kpc = 3.0857E16 km, 1 Gyr = 3.1557E16 sec
    # 1 kpc/Gyr = 3.0857/3.1557 = 0.977818
    kms2kpcGyr = 3.0857/3.1557
    kpcGyr2kms = 3.1557/3.0857
    npar = 8
    arr[:,npar-3:npar] = arr[:,npar-3:npar]*kpcGyr2kms


    # Remove the mean position and velocity
    for i in range(2,8):
        arr[:,i] = arr[:,i]-np.mean(arr[:,i])


    # ADDING INITIAL CONDITIONS
    # [mass,soft,x,y,z,Vx,Vy,Vz]
    arr[:,0] = 0.   # mass=0, test particles
    arr[:,1] = 0.   # no softening parameter
    arr[:,2] = arr[:,2] + x0
    arr[:,3] = arr[:,3] + y0
    arr[:,4] = arr[:,4] + z0
    arr[:,5] = arr[:,5] + vx0
    arr[:,6] = arr[:,6] + vy0
    arr[:,7] = arr[:,7] + vz0


    # FIRST PARTICLE IS THE GALAXY
    #arr(0,*) = [mass,rscale,0.,0.,0.,0.,0.,0.]
    arr = np.vstack(([mass,soft,x0,y0,z0,vx0,vy0,vz0], arr))

    # PRINTING PARAMETERS TO THE SCREEN

    # GETTING STRINGS
    if dhalo is None:
        dhalo=13.0
    strvxvyvz = str(vx0)+'  '+str(vy0)+'  '+str(vz0)
    strxyz = str(x0)+'  '+str(y0)+'  '+str(z0)
    strsoft = str(soft)
    strmmass = str(mass)

    # PRINTING INFO
    if noprint is None:
        # Printing the param inputs
        print('')
        print('PARAM INPUTS')
        print('-----------------------------------------------------')
        print('XYZ    = ',strxyz,'   kpc')
        print('VxVyVz = ',strvxvyvz,'   km/s')
        print('Soft   = ',strsoft,'   kpc')
        print('Mass   = ',strmmass,' Msun')
        print('-----------------------------------------------------')
    
    return arr

def rotate(coords,n1,n2,n3,t1):
    '''
    rotate(coords,n1,n2,n3,t1)
    
    Rotates position and velocity by angle t1 about a unit vector with
    coordinates (n1,n2,n3).
    
    Inputs
    ------
    coords:   2-D numpy array with [x,y,z,vx,vy,vz] for each body
    n1,n2,n3: floats or 1-D arrays of length nbodies, the x, y, and z
              coordinates of the vector to rotate about
    t1:       float or 1-D array of length nbodies, the rotation angle
    
    Returns
    -------
    coords1: 2-D numpy array with [x,y,z,vx,vy,vz] for each rotated body
    '''
    s1 = np.sin(t1)
    c1 = np.cos(t1)
    
    nbodies = len(coords[:,0])
    
    # 3-D rotation matrix
    R = np.zeros((nbodies,3,3))
    R[:,0,0] = c1+n1**2*(1-c1)
    R[:,0,1] = n1*n2*(1-c1)-n3*s1
    R[:,0,2] = n1*n3*(1-c1)+n2*s1
    R[:,1,0] = n1*n2*(1-c1)+n3*s1
    R[:,1,1] = c1+n2**2*(1-c1)
    R[:,1,2] = n2*n3*(1-c1)-n1*s1
    R[:,2,0] = n1*n3*(1-c1)-n2*s1
    R[:,2,1] = n2*n3*(1-c1)+n1*s1
    R[:,2,2] = c1+n3**2*(1-c1)

    # New array with rotated values
    coords1 = np.zeros(np.shape(coords))
    
    # x,y,z
    xyz0 = np.reshape(np.repeat(coords[:,:3],3),[nbodies,3,3])
    coords1[:,:3] = R[:,:,0]*xyz0[:,0,:]+R[:,:,1]*xyz0[:,1,:]+\
                    R[:,:,2]*xyz0[:,2,:]
    # vx,vy,vz
    xyz0 = np.reshape(np.repeat(coords[:,3:],3),[nbodies,3,3])
    coords1[:,3:] = R[:,:,0]*xyz0[:,0,:]+R[:,:,1]*xyz0[:,1,:]+\
                    R[:,:,2]*xyz0[:,2,:]
    
    return coords1

def todisk(coords):
    x = coords[:,0]
    y = coords[:,1]
    z = coords[:,2]
    vx = coords[:,3]
    vy = coords[:,4]
    vz = coords[:,5]
    
    Lx = y*vz-z*vy
    Ly = z*vx-x*vz
    Lz = x*vy-y*vx
    t0 = np.arccos(Lz/np.sqrt(Lx**2+Ly**2+Lz**2))
    
    ax = Lx*np.sin(t0/10)/np.sin(t0)
    ay = Ly*np.sin(t0/10)/np.sin(t0)
    az = Lz*np.cos(t0/10)/np.cos(t0)
    
    t1 = 9*t0/10
    
    bx = Ly*az-Lz*ay
    by = Lz*ax-Lx*az
    bz = Lx*ay-Ly*ax
    b = np.sqrt(bx**2+by**2+bz**2)
    
    n1 = bx/b
    n2 = by/b
    n3 = bz/b
    
    coords1 = rotate(coords,n1,n2,n3,t1)
    
    return coords1

def rotategal(coords,theta,phi):
    phi2 = phi+np.pi/2
    n1 = np.cos(phi2)
    n2 = np.sin(phi2)
    
    coords1 = rotate(coords,n1,n2,0,theta)
    
    return coords1

def makedisk(input0,nbody,noprint=None,stp=None,dhalo=None,theta=0,phi=0):


    #+
    # This program runs an N-test-body simulation of a
    # galaxy with a plummer distribution of stars that
    # orbits the Milky Way.
    #
    # INPUTS:
    #  input    Input parameters [mass,soft,x,y,z,vx,vy,vz] of the galaxy CM
    #  =nbody   Number of bodies to use in N-test-body simulation
    #  /noprint Don't print out the parameters
    #  /stp     Stop at the end of the program
    #
    # OUTPUTS:
    #  arr      The array of parameters for the galaxy.
    #
    # CALLS:
    #  plummer   To get the distribution of the N-bodies
    #
    # By David Nidever   March 2006
    #-

    # Not enough inputs
    if len(input0)==0:
        print('Syntax - makegal,input,arr,nbody=nbody,noprint=noprint,stp=stp')
        return

    # INITIAL CONDITIONS
    #input = [mass,soft,xstart,ystart,zstart,vxstart,vystart,vzstart]
    mass = input0[0]
    soft = input0[1]
    x0 = input0[2]
    y0 = input0[3]
    z0 = input0[4]
    vx0 = input0[5]
    vy0 = input0[6]
    vz0 = input0[7]


    # GETTING THE PLUMMER DISTRIBUTION OF STARS
    #mass = 1e7
    #rscale = 0.9*(mass/1d9)^(1./3.)

    print('USING {} STARS'.format(nbody))
    #print('Getting Plummer distribution of stars')
    arr = plummer(n=nbody,mass=mass,rscale=soft,soft=1)


    # VELOCITIES FROM PLUMMER ARE IN KPC/GYR
    # CONVERT TO KM/S FOR GALSAT INPUT

    # CONVERTING VELOCITIES (kpc/Gyr -> km/s)
    # 1 kpc = 3.0857E16 km, 1 Gyr = 3.1557E16 sec
    # 1 kpc/Gyr = 3.0857/3.1557 = 0.977818
    kms2kpcGyr = 3.0857/3.1557
    kpcGyr2kms = 3.1557/3.0857
    npar = 8
    arr[:,npar-3:npar] = arr[:,npar-3:npar]*kpcGyr2kms


    # Remove the mean position and velocity
    for i in range(2,8):
        arr[:,i] = arr[:,i]-np.mean(arr[:,i])

    arr[:,2:] = todisk(arr[:,2:])
    arr[:,2:] = rotategal(arr[:,2:],theta,phi)

    # ADDING INITIAL CONDITIONS
    # [mass,soft,x,y,z,Vx,Vy,Vz]
    arr[:,0] = 0.   # mass=0, test particles
    arr[:,1] = 0.   # no softening parameter
    arr[:,2] = arr[:,2] + x0
    arr[:,3] = arr[:,3] + y0
    arr[:,4] = arr[:,4] + z0
    arr[:,5] = arr[:,5] + vx0
    arr[:,6] = arr[:,6] + vy0
    arr[:,7] = arr[:,7] + vz0


    # FIRST PARTICLE IS THE GALAXY
    #arr(0,*) = [mass,rscale,0.,0.,0.,0.,0.,0.]
    arr = np.vstack(([mass,soft,x0,y0,z0,vx0,vy0,vz0], arr))

    # PRINTING PARAMETERS TO THE SCREEN

    # GETTING STRINGS
    if dhalo is None:
        dhalo=13.0
    strvxvyvz = str(vx0)+'  '+str(vy0)+'  '+str(vz0)
    strxyz = str(x0)+'  '+str(y0)+'  '+str(z0)
    strsoft = str(soft)
    strmmass = str(mass)

    # PRINTING INFO
    if noprint is None:
        # Printing the param inputs
        print('')
        print('PARAM INPUTS')
        print('-----------------------------------------------------')
        print('XYZ    = ',strxyz,'   kpc')
        print('VxVyVz = ',strvxvyvz,'   km/s')
        print('Soft   = ',strsoft,'   kpc')
        print('Mass   = ',strmmass,' Msun')
        print('-----------------------------------------------------')
    
    return arr