import corner
import emcee
import numpy as np

def InclineDiskDist(ra, dec, pa, io, d0):
    """
    Calculate the distance to a position in the plane of the LMC stellar disk

    Parameters
    ------------
    ra: right ascension of position
    dec: declination of position
    pa: line of nodes of the LMC stellar disk
    io: inclination of the LMC stellar disk
    d0: distance to the LMC's center of mass 

    Returns
    --------
    dist: distance to the position
    """
    ra1 = np.copy(ra)
    dec1 = np.copy(dec)
    ra1 = np.radians(ra1)
    dec1 = np.radians(dec1)
    alph0 = np.radians(82.25) #ra of LMC center
    delt0 = np.radians(-69.50) #dec of LMC center
    sd = np.sin(delt0)
    cd = np.cos(delt0)
    io = np.radians(io) #inclination #25.86
    pa = np.radians(pa + 90) #np.radians(149.37) #position angle of line-of-nodes
    cosr = cd*np.cos(dec1)*np.cos(ra1-alph0)+sd*np.sin(dec1)
    sinrcosph = -np.cos(dec1)*np.sin(ra1-alph0)
    sinrsinph = cd*np.sin(dec1)-sd*np.cos(dec1)*np.cos(ra1-alph0)
    dist = d0*np.cos(io)/(np.cos(io)*cosr-np.sin(io)*np.cos(pa)*sinrsinph+np.sin(io)*np.sin(pa)*sinrcosph)
    return dist


def lnL(theta, ra, dec, plx, plx_err):
    """
    Log Likelihood function

    Parameters
    ------------
    theta: line of nodes, inclination and central distance; the parameters to be fit
    ra: right ascension 
    dec: declination 
    plx: parallax
    plx_err: error in parallax 

    Returns
    --------
    lnl: log likelihood function
    """
    pa, io, d0 = theta
    model = InclineDiskDist(ra, dec, pa, io, d0)
    inv_sig2 = np.reciprocal(np.square(plx_err))
    lnl = -0.5 * np.sum(np.multiply(np.square((plx - model)),inv_sig2) - np.log(inv_sig2/(2*np.pi)))
    return lnl


def lnPrior(theta):
    """
    Prior distribution 

    Parameters
    ------------
    theta: line of nodes, inclination and central distance; the parameters to be fit

    Returns
    --------
    zero or infinity 
    """
    pa, io, d0 = theta
    if 0 < pa < 360 and 0 < io < 90.0 and 0 < d0 < np.inf: ## try narrower ranges
        return 0.0
    return -np.inf


def lnProb(theta, ra, dec, plx, plx_err):
    """
    Log of the posterior distribution 

    Parameters
    ------------
    theta: line of nodes, inclination and central distance; the parameters to be fit
    ra: right ascension 
    dec: declination 
    plx: parallax
    plx_err: error in parallax 
    
    Returns
    --------
    Log of the posterior distribution
    """
    lp = lnPrior(theta)
    if not np.isfinite(lp):
        return -np.inf
    return lp + lnL(theta, ra, dec, plx, plx_err)


# pa_y = 149.23 # line of nodes from Yumi 
# io_y = 25.86 # inclination from Yumi 
# d0_g = 49.59 # distance from Pietrzynski (2019)


# #To run the MCMC:
# guess = [pa_y, io_y, d0_g] #guess parameters

# pos = guess + np.random.randn(100, 3) ## do each variable seperately with 100 walkers
# nwalkers, ndim = pos.shape

# sampler = emcee.EnsembleSampler(nwalkers, ndim, lnProb, args=(ra, dec, plx, plx_err))
# sampler.run_mcmc(pos, 1000, progress=True); #1000 samples

# #To create corner plot:
# flat_samples = sampler.get_chain(discard=100, thin=1, flat=True)
# corner.corner(flat_samples, labels=labels, truths=[pa_y, io_y, d0_g]);

# #To print the results:
# vals = [np.median(flat_samples[:,0]), np.median(flat_samples[:,1]), np.median(flat_samples[:,2])]
# errs = [utils.mad(flat_samples[:,0]), utils.mad(flat_samples[:,1]), utils.mad(flat_samples[:,2])]

# print('line of nodes: {} +/- {}'.format(vals[0], errs[0]))
# print('inclination: {} +/- {}'.format(vals[1], errs[1]))
# print('distance: {} +/- {}'.format(vals[2], errs[2]))