import numpy as np
import time
from galsat import galsat
from outputread import outputread
from makegal import makegal
from plotnbody import plotnbody

def lsmc_orbit(input0=None,nbody=None,movie=None,df=None,R0=None,vcirc=None,\
               noplot=None,step=None,ostep=None,dstep=None,fstep=None,\
               dhalo=None,minstep=None,maxstep=None,tstart=None):

    #+
    # This program runs an N-test-body simulation of a
    # galaxy with a plummer distribution of stars that
    # orbits the Milky Way.
    #
    # INPUTS:
    #  input0    Input parameters [t,mass,soft,x,y,z,vx,vy,vz]
    #  =nbody   Number of bodies to use in N-test-body simulation
    #  /movie   Show the movie
    #  /df      Use dynamical friction
    #  =R0      Galactocentric distance of sun
    #  =Vcirc   Circular rotation velocity of the sun
    #  /noplot  Don't plot anything
    #
    # OUTPUTS:
    #  output   Output parameters 
    #
    # CALLS:
    #  plummer   To get the distribution of the N-bodies
    #  galsat    To run the N-test-body simulation
    #  plotnbody To plot the results or show the movie
    # 
    #
    # By David Nidever   March 2006
    #-

    if dhalo is None:
        dhalo = 13.0
    if vcirc is None:
        vcirc = 220.
    if R0 is None:
        R0 = 8.5
    if df is None:
        df = 1
    if nbody is None:
        nbody=1000

    #######################################################
    ## STEP 1:  The current parameters
    #######################################################

    if tstart is None:
        tstart = 1.0

    ####################################
    ## For the LMC
    ####################################
    # From Kallivayalil et al.2006
    # Vx = -86+/-12, Vy = -268+/-11, Vz = 252+/-16 km/s   in GSR
    #Vxlmc = -86.0
    #Vylmc = -268.0
    #vzlmc = 252.0
    vxlmc = -57
    vylmc = -226
    vzlmc = 221

    # Calculating the galactocentric xyz coordinates for the LMC
    # Taken from Connors et al. 2006
    #R0 = 8.5
    #distlmc = 49.43
    #gllmc = 280.46              #current l,b
    #gblmc = -32.89
    #lbd2xyz,gllmc,gblmc,distlmc,xlmc,ylmc,zlmc,R0=R0,/noprin
    xlmc = -0.8
    ylmc = -41.6
    zlmc = -27

    lmass = 1.8e11#2d10
    lsoft = 0.9*(lmass/1e9)**(1./3.)  # rscale
    #lsoft = 1.0

    inplmc = [lmass,lsoft,xlmc,ylmc,zlmc,vxlmc,vylmc,vzlmc]


    ####################################
    ## For the SMC
    ####################################
    # From Kallivayalil et al.2007
    # Vx = -87+/-48, Vy = -247+/-42, Vz = 149+/-37 km/s   in GSR
    #Vxsmc = -87.0
    #Vysmc = -247.0
    #Vzsmc = 149.0
    vxsmc = 19
    vysmc = -153
    vzsmc = 153

    # Calculating the galactocentric xyz coordinates for the LMC
    # Taken from Connors et al. 2006
    #R0 = 8.5
    #distsmc = 57.02
    #glsmc = 302.79             #current l,b
    #gbsmc = -44.30
    #lbd2xyz,glsmc,gbsmc,distsmc,xsmc,ysmc,zsmc,R0=R0,/noprin
    xsmc = 15.3
    ysmc = -36.9
    zsmc = -43.3

    smass = 2.1e10#1e9   # SMC
    ssoft = 0.9*(smass/1e9)**(1./3.)  # rscale

    inpsmc = [smass,ssoft,xsmc,ysmc,zsmc,vxsmc,vysmc,vzsmc]


    #######################################################
    ## STEP 2:  Run single particles backwards to get
    ##          initial condition
    #######################################################

    print('Running the two single particles backwards')

    # Making the intial array
    input1 = np.vstack((inplmc,inpsmc))

    # Running GALSAT
    galsat(input1,step=0.01,ostep=0.01,dstep=tstart,tmin=tstart,prog='galsat6',
           dhalo=dhalo,vcirc=vcirc,df=df,prefix='lsmc/run1')

    # Getting the initial conditions
#     nlast = n_elements(output1[*,0,0])-1
#     inplmc2 = reform(output1[nlast,0,1:8])
#     inpsmc2 = reform(output1[nlast,1,1:8])
    output1 = outputread('lsmc/run1.out')
    tlast = np.min(output1['TIME'])
    nlast = output1[(output1['TIME']==tlast)]
    inplmc2 = list(nlast[0])[2:]
    inpsmc2 = list(nlast[1])[2:]



    #######################################################
    ## STEP 3:  Make the galaxies
    #######################################################

    #nbody = 1000

    # Making the LMC parameters
    print('Making the LMC')
    lmcarr = makegal(inplmc2,nbody=nbody)

    # Making the SMC parameters
    print('Making the SMC')
    smcarr = makegal(inpsmc2,nbody=nbody)


    #######################################################
    ## STEP 4:  Run the full N-body simulation
    #######################################################


    # Combining the arrays
    input2 = np.vstack((lmcarr,smcarr))


    t0 = time.time()

    # SETTING GALSAT PARAMETERS
    if step is None:
        step = 0.001               # step size control parameter 
    if movie is not None:
        ostep = 0.01
    else:
        ostep = tstart  # output interval
    dstep = tstart                                           # diagnostic interval
    if fstep is None:
        fstep = 0.001             # fixed timestep

    # RUNNING GALSAT
    galsat(input2,step=step,ostep=ostep,dstep=dstep,tmin=tstart,df=df,\
           dhalo=dhalo,minstep=minstep,maxstep=maxstep,fstep=fstep,forward=1,\
           vcirc=vcirc,prog='galsat6',prefix='lsmc/run2')
    
    output2 = outputread('lsmc/run2.out')
    output3 = outputread('lsmc/run2.out',dt='arr')

    # Correct the time
    output2['TIME'] = output2['TIME']-tstart
    output3[:,:,0] = output3[:,:,0]-tstart

    print('{} seconds'.format(time.time()-t0))
    
    if movie is not None:
        plotnbody(arr=output3,notrail=True,afile='lsmc4000',dir0='lsmc4000')


#     # PLOTTING THE LAST SNAP (CURRENT)
#     nsnap = n_elements(output2[*,0,0])
#     if not keyword_set(noplot):plotnbody,arr=output2,/last,R0=R0

#     # SHOWING THE MOVIE
#     if keyword_set(movie):begin

#       ## Use the napo at the end
#       #colarr = reform(output2[nsnap-1,*,9])

#       ## Making reasonable colors, bound stars are white
#       #colarr2 = (colarr gt 0)*(colarr*30.+80) + (colarr eq 0)*255.
#       #colarr2(0) = 250   # satellite is red

#       colarr2=[fltarr(nbody)+150,fltarr(nbody)+250]
#       plotnbody,arr=output2,/movie,/notrail,ps=3,R0=R0,colarr=colarr2 #,/dotfirst,colarr=colarr2
    return output2