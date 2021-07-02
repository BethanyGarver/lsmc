from astropy.table import Table
import numpy as np

'''
outputread(filename,dt='Table')

Reads a galsat output file and converts the data into either an astropy Table or numpy array. Velocities are given in km/s

Inputs
------
filename: (str) The path of the file to read.
dt: (str, default 'Table') The type of output. Allows either 'Table' or 'arr'.

Returns
-------
t: if dt is 'Table', an astropy table with columns for 'TIME','INDEX','MASS','SOFT','X','Y','Z','VX','VY','VZ'
   if dt is 'arr', an array of dimensions (nSnaps,nBodies,9) with [t,mass,soft,x,y,z,vx,vy,vz] for each particle at each time
'''
def outputread(filename,dt='Table'):
    # Check if dt is an allowed value
    if dt!='Table' and dt!='arr':
        print("dt must be 'Table' or 'arr")
        return
    
    # Create Table
    if dt=='Table':
        t = Table(names=('TIME','INDEX','MASS','SOFT','X','Y','Z','VX','VY','VZ'))
    
    # Read file
    with open(filename) as f:
        lines = f.readlines()
        nBodies = int(lines[0])
        jump = nBodies+2
        nSnaps = int(len(lines)/jump)
        # Create array
        if dt=='arr':
            t = np.zeros((nSnaps,nBodies,9))
        
        # Run through snapshots
        for i in range(nSnaps):
            time = float(lines[jump*i+1])
            # Run through particles for each time
            for j in range(nBodies):
                index = int(j)
                line = lines[jump*i+j+2].split()
                for k in range(len(line)):
                    line[k] = float(line[k])
                if dt=='Table':
                    t.add_row([time,index]+line)
                else:
                    t[i,j,0] = time
                    t[i,j,1:9] = line
    f.close()
    
    # Convert velocities to km/s
    kpcGyr2kms = 3.1557/3.0857
    if dt=='Table':
        t['VX'] = t['VX']*kpcGyr2kms
        t['VY'] = t['VY']*kpcGyr2kms
        t['VZ'] = t['VZ']*kpcGyr2kms
    else:
        t[:,:,6:] = t[:,:,6:]*kpcGyr2kms
    return t