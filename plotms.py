import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm as LN
from outputread import outputread
from astropy import units as u
from astropy.coordinates import SkyCoord
import MagellanicStream as ms

def plotms(output):
    last = output[(output['TIME']==np.max(output['TIME']))]
    c_gc = SkyCoord(x=last['X'],y=last['Y'],z=last['Z'],unit=u.kpc,\
                    frame='galactocentric')
    c_ms = c_gc.transform_to(ms.MagellanicStream)
    ms_lg,ms_bg = c_ms.MSLongitude.degree,c_ms.MSLatitude.degree
    a = np.where(ms_lg>180)
    ms_lg[a] = ms_lg[a]-360
    plt.hist2d(ms_lg,ms_bg,norm=LN(),bins=100)
    plt.colorbar()
    plt.xlim(max(ms_lg),min(ms_lg))
    plt.xlabel('MS Longitude')
    plt.ylabel('MS Latitude')
    plt.show()