from makegal import makedisk,makegal
from galsat import galsat
from outputread import outputread
import numpy as np
from dlnpyutils.coords import xyz2lbd
from astropy import units as u
from astropy.coordinates import SkyCoord

def lmcinc(theta,phi,nbody=1000,tstart=1.0,lmc1=None,prefix=None,
           returnlmc1=False):
    '''
    Runs galsat for a disk with the LMC's current position with initial
    inclination defined by theta and phi
    
    Parameters
    ----------
    theta:      polar angle
    phi:        azimuthal angle
    nbody:      number of test particles in the disk (default 1000)
    tstart:     amount of time in Gyr to run the simulation back from the 
                present (default 1.0)
    lmc1:       input lmc coordinates at start time to skip the step where 
                galsat runs the single particle back in time
    prefix:     start of file name for input and output files
    returnlmc1: option to return lmc coordinates at start time in addition to 
                usual output
                
    Returns
    -------
    last: coordinates of particles at present time, including l, b, right 
          ascension, declination, and distance
    lmc1: (if returnlmc1 is True) coordinates of lmc center at start time
    '''
    
    # Run LMC from present back in time
    if lmc1 is None:
        lmass = 1.8e11#2d10
        lsoft = 0.9*(lmass/1e9)**(1./3.)  # rscale
        lmc0 = [lmass,lsoft,-0.8,-41.6,-27,-57,-226,221]
        galsat(lmc0,step=0.01,ostep=tstart,dstep=1.0,tmin=tstart,
               prog='galsat6',prefix='lmcinc/lmc0')
        a = outputread('lmcinc/lmc0.out',dt='arr')
        lmc1 = a[1,0,1:]
    
    # Make inclined disk
    arr0 = makedisk(lmc1,nbody,theta=theta,phi=phi)
    
    # Set file name for input and output files
    if prefix is None:
        prefix = 'lmcinc/lmc_t{:04}p{:04}n{:06}s{:03}'\
                 .format(round(theta*1000),round(phi*1000),nbody,
                         round(tstart*10))
    
    # Run galsat forward in time
    galsat(arr0,tmin=tstart,step=0.001,ostep=tstart,forward=True,
           prog='galsat6',prefix=prefix)
    
    # Read output file
    out0 = outputread(prefix+'.out')
    # Present time only
    last = out0[(out0['TIME']==max(out0['TIME']))]
    # l, b, and distance
    l,b,d = xyz2lbd(last['X'],last['Y'],last['Z'])
    last['L'] = l
    last['B'] = b
    last['DIST'] = d
    
    # Right ascension and declination
    c_gal = SkyCoord(l=l*u.degree,b=b*u.degree,frame='galactic')
    c_icrs = c_gal.icrs
    last['RA'] = c_icrs.ra
    last['DEC'] = c_icrs.dec
    
    if returnlmc1:
        return last,lmc1
    else:
        return last