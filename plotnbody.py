from outputread import outputread
import numpy as np
import matplotlib.pyplot as plt
import os

def plotnbody(file=None,over=None,twoplot=None,arr=None,stp=None,
              xr=None,yr=None,zr=None,rr=None,tr=None,iso=None,movie=None,
              notrail=None,anim=True,afile=None,dir0=None,colarr=None,
              first=None,ps=None,dotfirst=None,last=None,R0=None):

    #+
    #
    # PLOTNBODY
    #
    # This program plots the nbody orbits
    #
    # INPUTS:
    #  file      File containing the galsat output data
    #  dfile     File containing the galsat diagnostic data
    #  /over
    #  /twoplot  Only two panels
    #  arr=arr   Input the galsat output array
    #  /stp      Stop at end
    #  xr=xr     X-axis range
    #  yr=yr     Y-axis range
    #  zr=zr     Z-axis range
    #  rr=rr     Radius-axis range
    #  tr=tr     Time-axis range
    #  /iso      Make plots isotropic, aspect ratio=1
    #  /movie    Flip through the timesteps
    #  /notrail  If /movie set then don't show the trail for each body
    #  /anim     Make an animation, save snaps to files
    #  afile=afile Suffix for animation files
    #  colarr=colarr  Array of colors for the bodies
    #  /first   Only plot the first snap
    #  /last      Only plot the last snap
    #  /dotfirst  Make a dot for the first particle (the main body)
    #
    #  By David Nidever
    #-

    # Not enough inputs
    if file is None and arr is None:
        print('Syntax - plotnbody,file,dfile,over=over,twoplot=twoplot,arr=arr,')
        print('           stp=stp,xr=xr,yr=yr,zr=zr,rr=rr,tr=tr,iso=iso,')
        print('           movie=movie,notrail=notrail,anim=anim,afile=afile,')
        print('           dir=dir,colarr=colarr,first=first,ps=ps,')
        print('           dotfirst=dotfirst,last=last,R0=R0')
        return

    if anim is not None:
        movie=1

    # importing nbody data
    if arr is None:
        arr = outputread(file,dt='arr')

    if R0 is None:
        R0=8.5

    origarr = arr

    # if there are softening parameters then get rid of them.
    if len(arr[0,0,:])>8:
        dum = arr
        arr = arr[:,:,0:8]*0.
        arr[:,:,0:2] = dum[:,:,0:2]
        arr[:,:,2:8] = dum[:,:,3:9]
        dum = 0.

    ns = len(arr[:,0,0])
    nb = len(arr[0,:,0])
    # R for LMC
    r = np.zeros((ns,nb))
    for i in range(nb):
        for j in range(ns):
            r[j,i] = np.linalg.norm(arr[j,i,2:5])

    # Getting the index for t=0
    dum = np.abs(arr[:,0,0])
    dum = np.argmin(dum)

#     erase  # erase plot window

#     setdisp

    # colors for my screen
    red = 250
    lred = 210
    green = 190000
    orange = 310000
    yellow = 450000
    blue = -25000
    lblue = -15000
    purple = -20000
    white = -1

    # setting for postscript
    if False:#!d.name=='PS' or keyword_set(anim):
        loadct,39
        black=0
        purple=30
        blue=60
        aqua=80
        green=155   #135
        yellow=200
        orange=225
        white=0
        red=250
        lred=240
        #backg=white
        backg=yellow
        coarr = [green,aqua,yellow,blue,purple,lred,orange]

    thk = 2      # line thickness

    #linestyle: 0-solid line, 1-dotted, 2-dashed, 3-dot dash, 4-dot dot
    #          dot dash, 5-long dashes

#     psym8,0.7  #1.5

    co1 = green
    co2 = red
    co3 = orange
    co4 = blue

    # setting the ranges
    dx=np.max(arr[:,:,2])-np.min(arr[:,:,2])
    dy=np.max(arr[:,:,3])-np.min(arr[:,:,3])
    dz=np.max(arr[:,:,4])-np.min(arr[:,:,4])
    dr=np.max(r)-np.min(r)
    dt=np.max(arr[:,:,0])-np.min(arr[:,:,0])
    if xr is None:
        xr=[np.min(arr[:,:,2])-0.2*dx,np.max(arr[:,:,2])+0.2*dx]
    if yr is None:
        yr=[np.min(arr[:,:,3])-0.2*dy,np.max(arr[:,:,3])+0.2*dy]
    if zr is None:
        zr=[np.min(arr[:,:,4])-0.2*dz,np.max(arr[:,:,4])+0.2*dz]
    if rr is None:
        rr=[np.min(r)-0.2*dr,np.max(r)+0.2*dr]
    if tr is None:
        tr=[np.min(arr[:,:,0])-0.2*dt,np.max(arr[:,:,0])]

    if iso is not None:
        lo = np.min([xr[0],yr[0],zr[0]])
        hi = np.max([xr[1],yr[1],zr[1]])
        xr = [lo,hi]
        yr = [lo,hi]
        zr = [lo,hi]

    co1 = white
    thk = 0.5
    if ps is None:
        ps=8

#     if !d.name=='PS':
#         thk = 4.0

    dotsize = 1.5
    dotco = 250

    # MOVIE
    if movie is not None:
        filelist = []

        co2 = white
        if notrail is not None:
            co1=white
        else:
            co1=red

        nsnap = len(arr[:,0,0])
#         psym8,0.7
        if anim is not None:
            if dir0 is None:
                dir0='plotnbody'
            os.system('mkdir {}'.format(dir0))

        # Looping through snaps
        for j in range(nsnap):
            current = arr[j,0,0]
            
            plt.figure(figsize=(8.5,8))

            # Opening file or erasing
            if anim is not None:
                #if file_search(dir0,/test_dir)=='':file_mkdir,dir0
                if afile is None:
                    afile = 'plotnbody'
                num = '{:04}'.format(j)
                #num = strtrim(long(j),2)
                #if num lt 10 then num='0'+num
                #if num lt 100 then num='0'+num
                #if num lt 1000 then num='0'+num
                filename = dir0+'/'+afile+'_'+num+'.png'
                filelist.append(filename)
                #if j mod 25==0 then ps_open,filename,/color
#                 !p.font = 0
#                 loadct,39,/silent
#                 ps_open,filename,/color,/encap,thick=5
#                 device,/inches,xsize=8.5,ysize=8.0

            # Z vs X  (upper-left)
            #!p.multi=[4,2,2]
            #plot,dist(5),/nodata,xtit='XGC',ytit='ZGC',xr=xr,yr=zr,xs=1,ys=1
            plt.subplot(2,2,1)
            plt.plot(arr[j,:,2],arr[j,:,4],'r,')
            plt.xlabel('XGC')
            plt.ylabel('ZGC')
            plt.xlim(xr)
            plt.ylim(zr)
            if notrail is None:
                for i in range(nb):
                    plt.plot(arr[:j,i,2],arr[:j,i,4],'k-',linewidth=1)
#                 ps1 = ps & co1 = 0 & symsize = 1.0
#                 if dotfirst is not None and i==0:
#                     ps1=8
#                     co1=dotco
#                     symsize=dotsize

#                 # Use the input color
#                 if colarr is not None:
#                     co1 = colarr(i)

#                 #Overplotting the last position (to erase the position's dot)
#                 #oplot,[arr((j-1)>0,i,2)],[arr((j-1)>0,i,4)],linestyle=0,co=0,thick=thk,ps=ps1,symsize=symsize  

#                 #Overplotting the Current position 
#                 oplot,[arr(j,i,2)],[arr(j,i,4)],linestyle=0,co=co1,thick=thk,ps=ps1,symsize=symsize

#                 # Overplotting the trail
#                 if not keyword_set(notrail):$
#                     oplot,[arr(0:j,i,2)],[arr(0:j,i,4)],linestyle=0,co=co2,thick=thk

#                 # Overplotting the GC and Sun
#                 oplot,[0],[0],ps=1  #GC
#                 oplot,[-R0],[0],ps=2  #SUN
#             end
#             # Overplot the body last so you can see it
#             if keyword_set(dotfirst):begin
#                 oplot,[arr(j,0,2)],[arr(j,0,4)],linestyle=0,co=co1,thick=thk,co=dotco,ps=8,symsize=dotsize        # current position
#                 #if not keyword_set(notrail) then oplot,[arr(0:j,0,0)],[r(0:j,0)],linestyle=0,co=co2,thick=thk
#             endif

            # R vs T (upper-right)
            plt.subplot(2,2,2)
            plt.plot(arr[j,:,0],r[j,:],'r,')
            plt.xlabel('t')
            plt.ylabel('RGC')
            plt.xlim(tr)
            plt.ylim(rr)
            if notrail is None:
                for i in range(nb):
                    plt.plot(arr[j,:,0],r[j,:],'k-',linewidth=1)
#             !p.multi=[3,2,2]
#             plot,dist(5),/nodata,xtit='t',ytit='RGC',xr=tr,yr=rr,xs=1,ys=1
#             for i=0.,nb-1:begin
#                 ps1 = ps & co1 = 0 & symsize = 1.0
#                 if keyword_set(dotfirst) and i==0:begin
#                     ps1=8
#                     co1=dotco
#                     symsize=dotsize
#                 endif

#                 # Use the input color
#                 if keyword_set(colarr):co1 = colarr(i)

#                 # Overplotting the last position
#                 #oplot,[arr((j-1)>0,i,0)],[r((j-1)>0,i)],linestyle=0,co=0,thick=thk,ps=ps1,symsize=symsize

#                 # Overplotting the current position
#                 oplot,[arr(j,i,0)],[r(j,i)],linestyle=0,co=co1,thick=thk,co=co1,ps=ps1,symsize=symsize
#                 if not keyword_set(notrail):oplot,[arr(0:j,i,0)],[r(0:j,i)],linestyle=0,co=co2,thick=thk
#             end
#             # Overplot the body last so you can see it
#             if keyword_set(dotfirst):begin
#                 oplot,[arr(j,0,0)],[r(j,0)],linestyle=0,co=co1,thick=thk,co=dotco,ps=8,symsize=dotsize        # current position
#                 #if not keyword_set(notrail) then oplot,[arr(0:j,0,0)],[r(0:j,0)],linestyle=0,co=co2,thick=thk
#             endif

            # Y vs X (lower-left)
            plt.subplot(2,2,3)
            plt.plot(arr[j,:,2],arr[j,:,3],'r,')
            plt.xlabel('XGC')
            plt.ylabel('YGC')
            plt.xlim(xr)
            plt.ylim(yr)
            if notrail is None:
                for i in range(nb):
                    plt.plot(arr[:j,i,2],arr[:j,i,3],'k-',linewidth=1)
#             !p.multi=[2,2,2]
#             plot,dist(5),/nodata,xtit='XGC',ytit='YGC',xr=xr,yr=yr,xs=1,ys=1
#             for i=0.,nb-1:begin
#                 ps1 = ps & co1 = 0 & symsize = 1.0
#                 if keyword_set(dotfirst) and i==0:begin
#                     ps1=8
#                     co1=dotco
#                     symsize=dotsize
#                 endif

#                 # Use the input color
#                 if keyword_set(colarr):co1 = colarr(i)

#                 #Overplotting last position
#                 #oplot,[arr((j-1)>0,i,2)],[arr((j-1)>0,i,3)],linestyle=0,co=0,thick=thk,ps=ps1,symsize=symsize

#                 # Overplotting the current position
#                 oplot,[arr(j,i,2)],[arr(j,i,3)],linestyle=0,co=co1,thick=thk,ps=ps1,symsize=symsize

#                 if not keyword_set(notrail):oplot,[arr(0:j,i,2)],[arr(0:j,i,3)],linestyle=0,co=co2,thick=thk
#                 oplot,[0],[0],ps=1  # GC
#                 oplot,[-R0],[0],ps=2  #SUN
#             end
#             # Overplot the body last so you can see it
#             if keyword_set(dotfirst):begin
#                 oplot,[arr(j,0,2)],[arr(j,0,3)],linestyle=0,co=co1,thick=thk,co=dotco,ps=8,symsize=dotsize        # current position
#                 #if not keyword_set(notrail) then oplot,[arr(0:j,0,0)],[r(0:j,0)],linestyle=0,co=co2,thick=thk
#             endif

            # Z vs Y  (lower-right)
            plt.subplot(2,2,4)
            plt.plot(arr[j,:,3],arr[j,:,4],'r,')
            plt.xlabel('YGC')
            plt.ylabel('ZGC')
            plt.xlim(yr)
            plt.ylim(zr)
            if notrail is None:
                for i in range(nb):
                    plt.plot(arr[:j,i,3],arr[:j,i,4],'k-',linewidth=1)
#             !p.multi=[1,2,2]
#             plot,dist(5),/nodata,xtit='YGC',ytit='ZGC',xr=yr,yr=zr,xs=1,ys=1
#             for i=0.,nb-1:begin
#                 ps1 = ps & co1 = 0 & symsize = 1.0
#                 if keyword_set(dotfirst) and i==0:begin
#                     ps1=8
#                     co1=dotco
#                     symsize=dotsize
#                 endif

#                 # Use the input color
#                 if keyword_set(colarr):co1 = colarr(i)

#                 # Overplotting the last position
#                 #oplot,[arr((j-1)>0,i,3)],[arr((j-1)>0,i,4)],linestyle=0,co=0,thick=thk,ps=ps1,symsize=symsize

#                 # Overplotting the current position
#                 oplot,[arr(j,i,3)],[arr(j,i,4)],linestyle=0,co=co1,thick=thk,ps=ps1,symsize=symsize

#                 if not keyword_set(notrail):oplot,[arr(0:j,i,3)],[arr(0:j,i,4)],linestyle=0,co=co2,thick=thk
#                 oplot,[0],[0],ps=1  # GC
#                 oplot,[0],[0],ps=2  #SUN
#             end
#             # Overplot the body last so you can see it
#             if keyword_set(dotfirst):begin
#                 oplot,[arr(j,0,3)],[arr(j,0,4)],linestyle=0,co=co1,thick=thk,co=dotco,ps=8,symsize=dotsize        # current position
#                 #if not keyword_set(notrail) then oplot,[arr(0:j,0,0)],[r(0:j,0)],linestyle=0,co=co2,thick=thk
#             endif

#             # Closing file or waiting
#             if keyword_set(anim):begin
#                 if !d.name=='PS':ps_close
#                 ps2gif,filename+'.eps',/eps
#             endif else begin
#                 #if keyword_set(notrail) then wait,0.01 else wait,0.04
#             endelse

#             #stop

#         end # for j

            plt.suptitle('Time = {} Gyr'.format(round(current,3)),\
                         backgroundcolor='w')
            plt.savefig(filename)
            if j!=ns-1:
                plt.close()

        if anim is not None:
            if os.path.exists('imagelist.txt'):
                os.remove('imagelist.txt')
            with open('imagelist.txt','w') as file:
                for item in filelist:
                    file.write('%s\n' % item)
            os.system('convert @imagelist.txt {}.gif'.format(afile))

#     else:

#         # NORMAL PLOTTING

#         if keyword_set(last):ind=ns-1 else ind = 0

#         # Z vs X  (upper-left)
#         !p.multi=[4,2,2]
#         plot,dist(5),/nodata,xtit='XGC',ytit='ZGC',xr=xr,yr=zr,xs=1,ys=1
#         for i=0.,nb-1:begin
#             if keyword_set(colarr):co1=colarr(i)

#             # All snaps
#             if not keyword_set(first) and not keyword_set(last):$
#                 oplot,arr(*,i,2),arr(*,i,4),linestyle=0,co=co1,thick=thk

#             # Only first or last snap
#             if keyword_set(first) or keyword_set(last):$
#                 oplot,[arr(ind,i,2)],[arr(ind,i,4)],co=co1,thick=thk,ps=ps

#             if not keyword_set(last):$
#                 oplot,[arr(cur,i,2)],[arr(cur,i,4)],linestyle=0,co=co1,thick=thk,ps=ps    #current position 
#             oplot,[0],[0],ps=1  #GC
#             oplot,[-R0],[0],ps=2  #SUN

#         end

#         # R vs T (upper-right)
#         !p.multi=[3,2,2]
#         plot,dist(5),/nodata,xtit='t',ytit='RGC',xr=tr,yr=rr,xs=1,ys=1
#         for i=0.,nb-1:begin
#             if keyword_set(colarr):co1=colarr(i)

#             # All snaps
#             if not keyword_set(first) and not keyword_set(last):$
#                 oplot,arr(*,i,0),r(*,i),linestyle=0,co=co1,thick=thk

#             # Only first or last snap
#             if keyword_set(first) or keyword_set(last):$
#                 oplot,[arr(ind,i,0)],[r(ind,i)],co=co1,thick=thk,ps=ps

#             if not keyword_set(last):$
#                 oplot,[arr(cur,i,0)],[r(cur,i)],linestyle=0,co=co1,thick=thk,ps=ps        # current position
#         end

#         # Y vs X (lower-left)
#         !p.multi=[2,2,2]
#         plot,dist(5),/nodata,xtit='XGC',ytit='YGC',xr=xr,yr=yr,xs=1,ys=1
#         for i=0.,nb-1:begin
#             if keyword_set(colarr):co1=colarr(i)

#             # All snaps
#             if not keyword_set(first) and not keyword_set(last):$
#                 oplot,arr(*,i,2),arr(*,i,3),linestyle=0,co=co1,thick=thk

#             # Only first or last snap
#             if keyword_set(first) or keyword_set(last):$
#                 oplot,[arr(ind,i,2)],[arr(ind,i,3)],co=co1,thick=thk,ps=ps

#             if not keyword_set(last):$
#                 oplot,[arr(cur,i,2)],[arr(cur,i,3)],linestyle=0,co=co1,thick=thk,ps=ps    #current position
#             oplot,[0],[0],ps=1  # GC
#             oplot,[-R0],[0],ps=2  #SUN
#         end

#         # Z vs Y  (lower-right)
#         !p.multi=[1,2,2]
#         plot,dist(5),/nodata,xtit='YGC',ytit='ZGC',xr=yr,yr=zr,xs=1,ys=1
#         for i=0.,nb-1:begin
#             if keyword_set(colarr):co1=colarr(i)

#             # All snaps
#             if not keyword_set(first) and not keyword_set(last):$
#                 oplot,arr(*,i,3),arr(*,i,4),linestyle=0,co=co1,thick=thk

#             # Only first or last snap
#             if keyword_set(first) or keyword_set(last):$
#                 oplot,[arr(ind,i,3)],[arr(ind,i,4)],co=co1,thick=thk,ps=ps

#             if not keyword_set(last):$
#                 oplot,[arr(cur,i,3)],[arr(cur,i,4)],linestyle=0,co=co1,thick=thk,ps=ps    #current position

#             oplot,[0],[0],ps=1  # GC
#             oplot,[0],[0],ps=2  #SUN
#         end

#     endelse

#     # diagnostics
#     if keyword_set(dfile):begin
#         diag = importdiag(dfile)
#         etot = reform(diag(1,*))
#         ekin = reform(diag(2,*))
#         epot = reform(diag(3,*))
#         steps = reform(diag(4,*))
#     endif

#     if keyword_set(stp):stop

#     !p.multi = 0
#     arr = origarr

#     end
