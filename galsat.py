import numpy as np
import os
#+
#
# GALSAT
#
# This program runs the c++ galsat program.  It's probably best to
# call it from a batch file. 
# NOTE: All velocities should be input in KM/S.  They will also be
#       output this way.  The velocity internal to the C++ galsat
#       program is kpc/Gyr, but galsat.pro and importnbody.pro 
#       automatically makes the right conversions.
#
# INPUT
#   input   Array of inputs [Nbodies, [mass,soft,x,y,z,Vx,Vy,Vz]]
#   step    Step size control parameter (default: step = 0.01)
#   dstep   Diagnostic interval (default: dstep = tmin)
#   ostep   Output interval (default: ostep = 0.01)
#   minstep Minimum timestep (positive), when using test particles
#   maxstep Maximum timestep (positive), when using test particles
#   fstep   Fixed timestep (positive), when using test particles
#   tmin    Minimum time to integrate to, ending time value (default tmin = 4.0)
#   prog    Name of the program to use (default prog='galsat8')
#   vcirc   Rotation velocity of the Milky Way disk (default vcirc = 220.0)
#   dhalo   Milky Way halo softening parameter (default ds = 13.0)
#   vhalo   Milky Way logarithmic halo Vhalo parameter (default vhalo
#              = 121.0).  Only for version galsat8 and later
#   /forward  Integrate forward instead of backwards (default: forward=0)
#   /mw     Do NOT use the Milky Way potential (default: mw=0)
#   /df     USE dynamical friction (default df=0)
#   /xtra   Output extra debugging information
#   /noprint Don't print anything to the screen
#
# OUTPUT
#   output  Array of galsat output. [Nsnap,Nbodies,[t,mass,soft,x,y,z,Vx,Vy,Vz]]
#   =diag   The diagnostics.  [Nsnap,[time, Etot, Ekin, Epot, Nsteps]]
#
# Created by David Nidever July 2005
#-

def galsat(input0,step=None,tmin=None,dstep=None,ostep=None,\
           prog=None,vcirc=None,dhalo=None,vhalo=None,forward=False,\
           mw=True,df=True,minstep=None,maxstep=None,fstep=None,\
           noprint=None,xtra=None,diag=None,stp=None,nodelete=None,\
           prefix='tgal'):

    # Bad Input Parameters
    if len(input0)==0: #n_params()==0 or len(input0)==0:
        print('syntax - galsat,input,output,step=step,dstep=dstep,ostep=ostep,')
        print('                prog=prog,vcirc=vcirc,dhalo=dhalo,vhalo=vhalo,tmin=tmin,foward=forward')
        print('                mw=mw,df=df,minstep=minstep,maxstep=maxstep,fstep=fstep')
        print('                noprint=noprint,diag=diag')
        return

    # Where am I?
    #dir = userdir()
    #dir = '~/'
    #bindir = dir+'nbody/'
    bindir = '~/Documents/galsat/src/'
    #bindir = '/net/halo/dln5q/galsat/'
    #bindir = '/Volumes/data/net/halo/dln5q/galsat/'

    # 1 kpc = 3.0857E16 km, 1 Gyr = 3.1557E16 sec
    # 1 kpc/Gyr = 3.0857/3.1557 = 0.977818
    kms2kpcGyr = 3.0857/3.1557
    kpcGyr2kms = 3.1557/3.0857

    # galsat versions
    if prog is None:
        prog='galsat8'  # 'galsat6'
    try:
        vers = int(prog[6])
    except:
        if prog=='galsat':
            vers=1
        else:
            vers=0
    if mw and vers < 1:
        print('MW potential not available before galsat1')
    if df and vers < 6:
        print('/df not available before galsat6')
    if forward and vers < 4:
        print('/forward not available before galsat4')
    if vcirc is not None and vers < 2:
        print('=vcirc not available before galsat2')
    if dhalo is not None and vers < 3:
        print('=dhalo not available before galsat3')
    if vhalo is not None and vers < 8:
        print('=vhalo not available before galsat3')
    if minstep is not None and vers < 4:
        print('=minstep not available before galsat4')
    if maxstep is not None and vers < 4:
        print('=maxstep not available before galsat4')
    if fstep is not None and vers < 4:
        print('=fstep not available before galsat4')

    # Defaults
    if tmin is None:
        tmin=4.0
    if step is None:
        step=0.01
    if dstep is None:
        dstep=tmin
    if ostep is None:
        ostep=0.01
    if vcirc is None:
        vcirc=220.
    if dhalo is None:
        dhalo=13.0
    if vhalo is None:
        vhalo=121.0
#     if forward is None:
#         forward=0
#     if mw is None:
#         mw=0
#     if df is None:
#         df=0
    if xtra is None:
        xtra = 0

    # # of Bodies
    input0 = np.array(input0)
    sz = input0.shape
    if len(sz)==1:
        nbody = 1
    else:
        nbody = sz[0]

    # CONVERTING VELOCITIES (km/s -> kpc/Gyr)
    input2 = input0
    # More than one particle
    if nbody > 1:
        if sz[1]==7: # mass, x, y, z, vx, vy, vz
            input2[:,4] = input2[:,4]*kms2kpcGyr
            input2[:,5] = input2[:,5]*kms2kpcGyr
            input2[:,6] = input2[:,6]*kms2kpcGyr
        elif sz[1]==8: # mass, soft, x, y, z, vx, vy, vz
            input2[:,5] = input2[:,5]*kms2kpcGyr
            input2[:,6] = input2[:,6]*kms2kpcGyr
            input2[:,7] = input2[:,7]*kms2kpcGyr
        else:
            raise ValueError('Format not understood')
    # One particle
    else:
        if sz[0]==7: # mass, x, y, z, vx, vy, vz
            input2[4] = input2[4]*kms2kpcGyr
            input2[5] = input2[5]*kms2kpcGyr
            input2[6] = input2[6]*kms2kpcGyr
        elif sz[0]==8: # mass, soft, x, y, z, vx, vy, vz
            input2[5] = input2[5]*kms2kpcGyr
            input2[6] = input2[6]*kms2kpcGyr
            input2[7] = input2[7]*kms2kpcGyr
        else:
            raise ValueError('Format not understood')
    vcirc_kpcgyr = vcirc*kms2kpcGyr
    vhalo_kpcgyr = vhalo*kms2kpcGyr


    # Getting inputs ready
    if tmin < 0.:
        strtmin = str(round(-tmin,2))
    else:
        strtmin = str(tmin).strip()
    strstep = str(step).strip()
    strdstep = str(dstep).strip()
    strostep = str(ostep).strip()
    #strvcirc = strtrim(vcirc,2)
    strvcirc = str(vcirc_kpcgyr).strip()
    strdhalo = str(dhalo).strip()
    #strvhalo = strtrim(vhalo,2)
    strvhalo = str(vhalo_kpcgyr).strip()
    strnbody = str(nbody).strip()

    # Setting flags
    flags = ''
    if forward:
        flags = flags+' -f'
    if not mw and vers >= 1:
        flags = flags+' -m'
    if df:
        flags = flags+' -y'
    if xtra!=0:
        flags = flags+' -x'

    # CREATING INPUT FILE
    #tfile = maketemp('tgal')
    filename = prefix+'.in'
    tfile = open(filename,'w')
    tfile.write(strnbody+'\n')
    tfile.write('0\n')
    if vers >= 2:
        tfile.write(strvcirc+'\n')    # Vcirc
    if vers >= 3:
        tfile.write(strdhalo+'\n')    # dhalo
    if vers >= 8:
        tfile.write(strvhalo+'\n')
    if nbody > 1:
        for i in range(nbody):
            newline = ''
            for j in input2[i,:]:
                newline = newline+str(j)+' '
            tfile.write(newline+'\n')
    #    for i=0.,nbody-1 do printf,unit,strtrim(input2(i,*),2)
    elif nbody==1:
        for j in input2:
            tfile.write(str(j)+' ')
        tfile
    #printf,unit,'2.0e10 ',strsoft,' ',strxyz,' ',strvxvyvz
    #printf,unit,'7.5e8 0.0 16.1556 2.27043 -5.88738  237.832 -42.3568 221.998'
    tfile.close()

    # Printing info
    if noprint is None:

        # Printing the param inputs
        print('PARAM INPUTS')
        print('-----------------------------------------------------')
        print('Prog = ',prog)
        print('Nbodies = ',strnbody)
        print('Tmin   = ',strtmin,'   Gyrs')
        if vers > 2:
            print('Vcirc  = ',str(vcirc).strip(),'   km/s')
        if vers > 3:
            print('Dhalo  = ',strdhalo,'   kpc')
        if vers >= 8:
            print('Vhalo  = ',str(vhalo).strip(),'   km/s')
        print('Step   = ',strstep)
        print('Ostep  = ',strostep,' Gyrs')
        print('Dstep  = ',strdstep,' Gyrs')
        if forward:
            print('Integrating FORWARDS')
        if not mw:
            print('NOT Using Milky Way Potential')
        else:
            print('Using Milky Way Potential')
        if df:
            print('Using Dynamical Friction')
        print('-----------------------------------------------------')

    # RUNNING GALSAT
    print(' ')
    #print('Running GALSAT program ... 3 sec.'
    print('Running GALSAT program ...')
    #spawn,'\rm '+tfile+'.out'
    cmd = '( '+bindir+prog+' -d '+strstep+' -o '+strostep+' -t '+strtmin+' -i -e '+strdstep
    if minstep is not None:
        cmd = cmd+' -l '+str(minstep).strip()
    if maxstep is not None:
        cmd = cmd+' -c '+str(maxstep).strip()
    if fstep is not None:
        cmd = cmd+' -s '+str(fstep).strip()
    cmd = cmd + flags+' < '+prefix+'.in > '+prefix+'.out )'# > & '+tfile+'.diag'

    #print,cmd

    #spawn,cmd,dum
    os.system(cmd)

#     #spawn,'./'+prog+' -d 0.01 -o 0.01 -t '+strtmin+' -i -e '+strtmin+' < '+tfile+'.in > '+tfile+'.out',dum

#     # GETTING THE OUTPUT
#     fileinfo = file_info(tfile+'.out')
#     fsize = fileinfo.size

#     # Everything in one file
#     if (fsize > 0):begin
#       output = importnbody(tfile+'.out')

#     # In separate files
#     endif else begin
#       galfiles = file_search('gal*.out')
#       nfiles = len(galfiles)

#       d = importnbody(galfiles(0))
#       sz = size(d)
#       output = dblarr(nfiles,sz(2),sz(3))

#       # Restoring the files
#       for i=0L,nfiles-1 do begin
#         d = importnbody(galfiles(i))
#         output(i,*,*) = d
#       end

#       # Sort by time
#       si = sort(output(*,0,0))
#       output = output(si,*,*)

#       stop

#     endelse

#     # Read in the diagnostic information
#     #  time  Etot  Ekin  Epot  Nsteps
#     readline,tfile+'.diag',diaglines,comment='#',count=ndiaglines
#     if keyword_set(xtra):begin
#       diaglines0 = diaglines
#       bd = where(stregex(diaglines,'internal data for',/boolean)==1,nbd)
#       if nbd > 0:remove,bd,diaglines
#       brk = where(stregex(diaglines,'for debugging',/boolean),nbrk)
#       lo = brk-1
#       hi = [brk[1:*]-2,len(diaglines)-1]
#       diag = fltarr(nbrk,nbody,18)
#       # loop through outputs
#       for i=0,nbrk-1 do begin
#         temp = diaglines[lo[i]:hi[i]]
#         sumline = temp[0]  # t, etot, ekin, epot, nsteps
#         sumlinearr = strsplit(sumline,' ',/extract)
#         time = float(sumlinearr[0])
#         temp = temp[2:*]  # cut out header lines
#         # mass, soft, x, y, z, vx, vy, vz, accx, accy, accz, jerkx, jerky, jerkz
#         #   epot, ekin, etot
#         tarr = strsplitter(temp,' ',/extract)
#         tarr = float(transpose(tarr))
#         diag[i,*,0] = time
#         diag[i,*,1:*] = tarr
#       endfor
#     endif

#     if keyword_set(stp):stop

#     # REMOVING TEMPORARY FILES
#     if not keyword_set(nodelete):file_delete,tfile+['.in','.out','.diag'],/allow

#     end

