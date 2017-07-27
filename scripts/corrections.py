#!/Users/fivos/pyvos/bin/python

import numpy as np
import sys


def pol_corr(center, d, pxd, polfac, nx,ny):
    """
    Polarization correction for all pixels ("pol" in Skinner et al. 2012):
    Correction for intensity modulation due to polarization of the beam
    needs to be done before angular integration.
	
	PARAMETERS
	--
	d	    : sample - detector distance [m]
	pxd		: pixel size [m]
	polfac	: polarization factor (percent of horizontal polarization)
	nx,ny	: number of pizels on the detector along x and y

	REFERENCES
	--
	L.B. Skinner et al. Nucl. Instr. 662, 61-70 (2012)
    """
    corrdat = np.zeros((nx,ny))#data
    xcnt, ycnt = center[0], center[1]
    for i in xrange(0, nx):
        inew = i - xcnt
        for j in xrange(0, ny):
            jnew = j - ycnt
            two_theta = np.arctan(np.sqrt(inew*inew + jnew*jnew) * pxd / d)
            if (inew == 0):
                phi = 0.0
            else:
                phi = np.arctan2(jnew, inew)
            alpha = polfac * (1. - np.sin(phi)**2. * np.sin(two_theta)**2.) + (1. - polfac) * (1. - np.cos(phi)**2. * np.sin(two_theta)**2.)
            corrdat[i, j]   = alpha#data[i, j] / alpha
    return corrdat

def SolAng_corr(q, pxd,  wlen, d):
    """
    This correction accounts for the solid angle covered by pixels at different theta
    ("geo" in Skinner et al., 2012)

    PARAMETERS
    --
	q		: array of the q values [A-1]
    d	    : sample - detector distance [m]
    pxd		: pixel size [m]
    nx,ny	: number of pizels on the detector along x and y
    wlen    : wavelength x-ray [A]

    REFERENCES
    --
	L.B. Skinner et al. Nucl. Instr. 662, 61-70 (2012)
    """
    Nbin = len(q)
    const = 1.#(pxd / d)**2.0
    SolAngCorr = np.zeros(Nbin)
    for i in xrange(0, Nbin):
        two_theta = q2th(q[i], wlen)
        SolAngCorr[i] = np.power(np.cos(two_theta), 3)
    return const*SolAngCorr

def q2th(q, wlen):
    """
    Convert momentum transfer to two-theta
    """
    th = np.arcsin(q*wlen / (4.0 * np.pi) )
    return th

