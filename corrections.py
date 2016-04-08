# ------------------------------------------------------------------------------
#	Python routines for the analysis of TIFF detector images at APS
#	including background subtraction, polarization and solid angle correction,
#	Compton scattering subtraction and normalization.
#
#	Daniel Schlesinger, Stockholm 2016-03-07
#   daniel.schlesinger@fysik.su.se
#
#	Last update: 2016-04-08
#	TO DO
#	- implement ring-finder, center finder
#	- masking
#
# ----------------------------------- Imports ----------------------------------
import sys
import numpy as np
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit
from scipy import signal
import csv
from PIL import Image
import scipy as sc
import gdal
from RingData import RingFit

# --------------------------------- Main program -------------------------------
def main():
	
	# -------------------------- Constants (CHECK!)---------------
	d			= 0.348			# sample - detector distance [m]
	wtodet		= 0.1			# kapton window to detector dist [m]
	pxd			= 200.0e-6		# pixel size [m]
	maxpx		= 2048			# pixel number (assuming square detector)
	wlen		= 0.123628      # wavelength x-ray [A]
	wlen_f		= 0.5			# fluorescence wavelength
	polfac		= 0.95			# polarzation factor
	qmax		= 35			# maximum q value [A^-1]
	dq			= 0.01			# q-step [A^-1]
	Ts			= 1.0			# q=0 transmission of sample
	Tc			= 1.0			# q=0 transmission of container
	mut_sam		= 6.0e-6		# droplet of 20mu
	dair		= 30			# air gap between sample and det.
	dkapton		= 2.0e-3		# kapton windon thickness
	#-------------------------------------------------------------
	
	# -------- read detector data --------
	datafile	= sys.argv[1]	# tif file with data image
	bkgrfile	= sys.argv[2]	# tif file with background image
	data		= read_TIFF(datafile)
	bkgr		= read_TIFF(bkgrfile)
	
	# --------- find symmetry axis > Will be replaced by the ring finder -------
	center		= np.zeros(2)	# image cnter (x0, y0) to be found
	#beta		= np.array([1031.3065, 1036.6479, 53.8183])   # ringfinder guess xcenter, ycenter, radius
	beta		= np.array([1031.3065, 1036.6479, 53.8183])
	
	# ------------- q-I(q) ---------------
	plot_detimg(data)
	r = RingFit(data)
	center[0], center[1], radius = r.fit_circle(beta)
	print "   xcenter, ycenter, radius: ", center[0], center[1], radius
	#sys.exit()
	datahisto = ang_int(data, center, pxd, maxpx, qmax, wlen, d, dq)
	#bkgrhisto = ang_int(bkgr, center, pxd, maxpx, qmax, wlen, d, dq)
	mplplot_2D(datahisto, qmax, dq)
	#plot_2D(bkgrhisto, qmax, dq)
	
	#generate auxiliary files:
	#F	= Fcorr(center, wlen_f, wlen, qmax, d, dq, pxd, maxpx, polfac, dair, dkapton, mut_sam)
	#write_Iq(F, 'Fcorr.dat', qmax, dq)
	#S	= Satt(mut_sam, qmax, wlen, d, dq)
	#aux	= np.ones(len(S))
	#write_Iq(aux / S, 'Satt.dat', qmax, dq)
	#C = Ccorr(wlen, qmax, dq, dair, dkapton)
	#KN = KleinNishina(wlen, qmax, dq)
	#plot_detimg(C)
	#write_Iq(C, 'C.dat', qmax, dq)
	
	
	# ---------- do corrections ----------
	#Is					= data - bkgr * Ts * Tc
	#Is_int				= ang_int(Is, center, pxd, maxpx, qmax, wlen, d, dq)
	#write_Iq(Is_int, 'normtest.dat', qmax, dq)
	#sys.exit()
	#Is_pc				= pol_corr(Is, center, d, pxd, polfac, maxpx)
	#Is_pc_int			= ang_int(Is_pc, center, pxd, maxpx, qmax, wlen, d, dq)
	#S					= Satt(mut_sam, qmax, wlen, d, dq)
	#geo					= SolAng_corr(data, pxd, qmax, wlen, d, dq)
	#Mf					= M(mut_det(wlen), mut_fil(wlen, dair, dkapton), wlen, qmax, dq)
	#Ix1					= Is_pc_int * geo / (Mf * S)
	#write_Iq(Ix1, 'firstterm-int2.dat', qmax, dq)
	#mplplot_2D(Ix1, qmax, dq)


# ---------------------------- Function definitions ----------------------------
# Skinner et al., 2012 = Skinner et al., Nuc. Inst. Meth. Phys. Res. A, 2012

def read_TIFF(filename):
	"""
	Read TIFF image and extract 2D detector image
	and return 2D numpy array (x,y)
	"""
	try:
		tif = gdal.Open(filename)
		tifArray = tif.ReadAsArray()
	except:
		print 'The file does not exist.'
		sys.exit(0)
	#print tifArray.shape, tifArray.size
	#plot_detimg(tifArray)
	return tifArray

def plot_detimg(tifArray):
	"""
	plot the detector image or array of same shape
	"""
	grid = tifArray
	plt.imshow(grid, interpolation=None)
	plt.show()
	return 0

def write_Iq(Iq, filename, qmax, dq):
	"""
	This function writes q-I(q) data to txt file
	"""
	f = open(filename, 'w')
	Nbin = int(qmax / dq)
	for i in xrange(1, Nbin):
		f.write(str(i*dq)+ "   "+str(Iq[i])+"\n")
	f.close
	return 0

def mask_badpxl():
	"""
	Bad pixel masking
	"""
	#calculate average background, then mask those pixels with high variance
	return 0

def pol_corr(data, center, d, pxd, polfac, maxpx):
	"""
	Polarization correction for all pixels ("pol" in Skinner et al. 2012):
	Correction for intensity modulation due to polarization of the beam
	needs to be done before angular integration.
	"""
	corrdat = data
	xcnt, ycnt = center[0], center[1]
	for i in xrange(0, maxpx):
		inew = i - xcnt
		for j in xrange(0, maxpx):
			jnew = j - ycnt
			two_theta = np.arctan(np.sqrt(inew*inew + jnew*jnew) * pxd / d)
			if (inew == 0):
				phi = 0.0
			else:
				phi	= np.arctan2(jnew, inew)
			#alpha	= (1.0 + np.cos(two_theta)*np.cos(two_theta) - polfac*np.cos(phi)*np.sin(two_theta)*np.sin(two_theta))/2.0
			alpha	= polfac * (1 - np.sin(phi)**2 * np.sin(two_theta)**2) + (1 - polfac) * (1 - np.cos(phi)**2 * np.sin(two_theta)**2)
			corrdat[i, j]	= data[i, j] / alpha
	return corrdat

def Satt(mut_sam, qmax, wlen, d, dq):
	"""
	Slab attenuation (Satt) as function of q, Skinner et al. 2012                 >TEST THIS A LITTLE MORE, STRANGE BEHAVIOR<
	"""
	Nbin = int(qmax / dq)
	Satt = np.zeros(Nbin)
	for i in xrange(0, Nbin):
		qval = i * dq
		two_theta = q2th(qval, wlen)
		ang_rad = np.pi*(two_theta / 180.0 + 0.5)
		denom	= np.sin(ang_rad)
		if (denom != 0):
			var = mut_sam / denom
		else:
			var = 0
		Satt[i] = np.exp(-var)*(np.exp(var - mut_sam) - 1) / (var - mut_sam)
	return Satt

def SolAng_corr(data, pxd, qmax, wlen, d, dq):
	"""
	This correction accounts for the solid angle covered by pixels at different theta
	("geo" in Skinner et al., 2012)
	"""
	Nbin = int(qmax / dq)
	const = (pxd / d)**2.0
	SolAngCorr = np.zeros(Nbin)
	for i in xrange(0, Nbin):
		qval = i * dq
		two_theta = q2th(qval, wlen)
		SolAngCorr[i] = np.power(np.cos(two_theta), 3)
	return const*SolAngCorr

def mut_fil(wlen, dair, dkapton):
	"""
	Reciprocal attenuation length times thickness mut_fil
	for filter. Include kapton & air.
	"""
	mu_air		= 0.00011486 + 0.00055966*wlen	# [cm^-1]; fit valid in range 70-100 keV
	mu_kapton	= 0.14901 + 0.55834*wlen		# [cm^-1]; fit valid in range 70-100 keV
	mut_fil		= dair*mu_air + dkapton*mu_kapton
	return mut_fil

def mut_det(wlen):
	"""
	Reciprocal attenuation length times thickness mut_det
	for detector.													>>THIS IS CURRENTLY UNKNOWN, CHECK!!!<<
	"""
	ddet		= 0.3						# cm
	mu_det		= 0.14901 + 0.55834*wlen
	mut_det		= ddet*mu_det
	return mut_det

def	M(wlen, qmax, dq, dair, dkapton):
	"""
	Attenuation corrections:
	Detector attenuation (mut_det)
	Filter attenuation   (mut_fil, includes air attenuation)
	Skinner et al., 2012
	"""
	Nbin	= int(qmax / dq)
	M		= np.zeros(Nbin)
	#det		= np.zeros(Nbin)
	#filter	= np.zeros(Nbin)
	for i in xrange(0, Nbin):
		qval = i * dq
		two_theta = q2th(qval, wlen)
		filter = np.exp(-mut_fil(wlen, dair, dkapton)/np.cos(two_theta)) / np.exp(-mut_fil(wlen, dair, dkapton))
		det	  = (1 - np.exp(-mut_det(wlen)/np.cos(two_theta))) / (1 - np.exp(-mut_det(wlen)))
	return det*filter

#------------ Compton correction

def mut_fil_array(wlen, dair, dkapton, qmax, dq):
	"""
	Reciprocal attenuation length times thickness mut_fil
	for filter. Include kapton & air. 
	Version with angle-dependent lambda
	"""
	Nbin		= int(qmax/dq)
	mut_fil		= np.zeros(Nbin)
	for i in xrange(0, Nbin):
		qval		= i * dq
		sin_two_th	= 4.0*np.pi*qval*wlen
		wlen_cs		= wlen + 0.048527 * sin_two_th**2.0
		mu_air		= 0.00011486 + 0.00055966*wlen_cs	# [cm^-1]; fit valid in range 70-100 keV
		mu_kapton	= 0.14901 + 0.55834*wlen_cs			# [cm^-1]; fit valid in range 70-100 keV
		mut_fil[i]	= dair*mu_air + dkapton*mu_kapton
	
	return mut_fil

def mut_det_array(wlen, qmax, dq):
	"""
	Reciprocal attenuation length times thickness mut_det
	for detector.													>>THIS IS CURRENTLY UNKNOWN, CHECK!!!<<
	Version with angle-dependent lambda
	"""
	Nbin		= int(qmax/dq)
	mut_det		= np.zeros(Nbin)
	ddet		= 0.3									# cm
	for i in xrange(0, Nbin):
		qval		= i * dq
		sin_two_th	= 4.0*np.pi*qval*wlen
		wlen_cs		= wlen + 0.048527 * sin_two_th**2.0
		mu_det		= 0.14901 + 0.55834*wlen_cs
		mut_det[i]	= ddet*mu_det
	
	return mut_det

def	M_array(wlen, qmax, dq, dair, dkapton):
	"""
	Attenuation corrections:
	Detector attenuation (mut_det)
	Filter attenuation   (mut_fil, includes air attenuation)
	Skinner et al., 2012
	"""
	Nbin	= int(qmax / dq)
	M		= np.zeros(Nbin)
	det		= np.zeros(Nbin)
	filter	= np.zeros(Nbin)
	#mf		= mut_fil_array(wlen, dair, dkapton, qmax, dq)
	#md		= mut_det_array(wlen, qmax, dq)
	for i in xrange(0, Nbin):
		qval = i * dq
		two_theta = q2th(qval, wlen)
		sin_two_th	= 4.0*np.pi*qval*wlen
		wlen_cs		= wlen + 0.048527 * sin_two_th**2.0
		#filter[i] = np.exp(-mut_fil_array(wlen, dair, dkapton, qmax, dq)/np.cos(two_theta)) / np.exp(-mut_fil_array(wlen, dair, dkapton, qmax, dq))
		#det[i]	  = (1 - np.exp(-mut_det_array(wlen, qmax, dq)/np.cos(two_theta))) / (1 - np.exp(-mut_det_array(wlen, qmax, dq)))
		filter[i] = np.exp(-mut_fil(wlen_cs, dair, dkapton)/np.cos(two_theta)) / np.exp(-mut_fil(wlen_cs, dair, dkapton))
		det[i]	  = (1 - np.exp(-mut_det(wlen_cs)/np.cos(two_theta))) / (1 - np.exp(-mut_det(wlen_cs)))
	
	return det*filter

def KleinNishina(wlen, qmax, dq):
	"""
	Klein-Nishina formula, eq. (12) in Skinner et al., 2012
	this goes into Compton scattering correction
	"""
	Nbin	= int(qmax / dq)
	h		= 6.62607004e-34  # m^2 kg / s
	c		= 2.99792458e+8   # m/s
	me		= 9.10938356e-31  # kg
	const	= h / (wlen * 1.0e-10 * me * c)						# OBS: wlen in Angstrom!
	KN		= np.zeros(Nbin)
	for i in xrange(0, Nbin):
		qval = i * dq
		two_theta = q2th(qval, wlen)
		P		= 1.0 / (1.0 + const * (1-np.cos(two_theta)))		# eq. (13)
		KN[i]	= (np.power(P,3) - np.power(P*np.sin(two_theta),2) + P)/(1+np.cos(two_theta)) # eq. (12)
	return KN

def Ccorr(wlen, qmax, dq, dair, dkapton):
	"""
	Compton scattering correction (inelastic scattering)
	eq. (12) in Skinner et al., 2012
	uses functions: mut_fil, mut_det, M, KleinNishina
	"""
	Nbin		= int(qmax / dq)
	M_CS		= M_array(wlen, qmax, dq, dair, dkapton)
	M_0			= M(wlen, qmax, dq, dair, dkapton)
	KN			= KleinNishina(wlen, qmax, dq)
	Ccorr		= KN * M_CS / M_0
	return Ccorr

#---------- end Compton.


def ang_int_simple(data, center, pxd, maxpx, qmax, wlen, d, dq):
	"""
	Angular integration of corrected detector image (binning into q-bins)
	"""
	qu		= pxd / d					# frequently used ratio
	two_pi	= 2.0 * np.pi				# frequently ...
	const	= 4.0 * np.pi / (wlen * dq)
	Nbin	= int(qmax / dq)			# number of q-bins
	histo	= np.zeros(Nbin)			# initialize q-bin histogram
	xcnt, ycnt = center[0], center[1]
	
	# loop over all pixels of the detector image
	for i in xrange(0, maxpx):
		inew = i - xcnt
		inew2 = inew*inew
		for j in xrange(0, maxpx):
			jnew = j - ycnt
			theta	= 0.5 * np.arctan(np.sqrt(inew2 + jnew*jnew) * qu)	# calculate angle
			qval = const * np.sin(theta)								# calc. q
			qbin = int(qval)											# binning acc. to q
			if ((qbin < Nbin) & (qbin > 0)):
				histo[qbin] += data[i, j] / (two_pi*qval)
	return histo

def ang_int(data, center, pxd, maxpx, qmax, wlen, d, dq):
	"""
	Angular integration of corrected detector image (binning into q-bins)
	"""
	qu		= pxd / d					# frequently used ratio
	two_pi	= 2.0 * np.pi				# frequently ...
	const	= 4.0 * np.pi / (wlen * dq)
	Nbin	= int(qmax / dq)			# number of q-bins
	histo	= np.zeros(Nbin)			# initialize q-bin histogram
	norm	= np.ones(Nbin)
	xcnt, ycnt = center[0], center[1]
	
	# loop over all pixels of the detector image
	for i in xrange(0, maxpx):
		inew = i - xcnt
		inew2 = inew*inew
		for j in xrange(0, maxpx):
			jnew = j - ycnt
			theta	= 0.5 * np.arctan(np.sqrt(inew2 + jnew*jnew) * qu)	# calculate angle
			qval = const * np.sin(theta)								# calc. q
			qbin = int(qval)											# binning acc. to q
			if ((qbin < Nbin) & (qbin > 0)):
				histo[qbin] += data[i, j]
				norm[qbin]	+= 1
					
	return histo / norm													# this won't be stable if dq too small.

def mplplot_2D(histo, qmax, dq):
	"""
	Line-plots of q-I(q) plots
	"""
	x = np.arange(0, qmax, dq)
	y = histo
	plt.plot(x, y)
	plt.show()
	return 0

def	q2th(q, wlen):
	"""
	Convert momentum transfer to theta
	"""
	th = np.arcsin(q*wlen / (4.0 * np.pi) )
	return th

def th2q(th, wlen):
	"""
	Convert theta to momentum transfer
	"""
	q = 4.0 * np.pi * np.sin(th) / wlen
	return q

def Fcorr(center, wlen_f, wlen_0, qmax, d, dq, pxd, maxpx, polfac, dair, dkapton, mut_sam):
	"""
	Fluorescence correction, q-dependent
	fluorescence wavelength wlen_f typically k-alpha						> OBS.: FUNCTION HAS NOT BEEN TESTED!!<
	"""
	ons   	= np.ones((maxpx, maxpx))
	polcorr	= pol_corr(ons, center, d, pxd, polfac, maxpx)
	pol		= ang_int_simple(polcorr, center, pxd, maxpx, qmax, wlen_0, d, dq)
	mff		= mut_fil(wlen_f, dair, dkapton, qmax, dq)
	mdf		= mut_det(wlen_f, qmax, dq)
	Mf		= M(mdf, mff, wlen_f, qmax, dq)
	mf0		= mut_fil(wlen_0, dair, dkapton, qmax, dq)
	md0		= mut_det(wlen_0, qmax, dq)
	M0		= M(md0, mf0, wlen_0, qmax, dq)
	Sf		= Satt(mut_sam, qmax, wlen_f, d, dq)
	S0		= Satt(mut_sam, qmax, wlen_0, d, dq)
	
	Fcorr	= Sf * Mf / (S0 * M0 * pol)
	
	return Fcorr


main()