import matplotlib.pyplot as plt
import numpy as np
import scipy as sp
from astropy.io import fits
from astropy.stats import sigma_clipped_stats
from photutils import daofind
from photutils import CircularAperture
from astropy.visualization import SqrtStretch
from astropy.visualization.mpl_normalize import ImageNormalize
import robust as rb


path = "/Users/nickedwards/Desktop/"

filename = path + "stars.fits"

# some standard variables
siglo=3
sighi=7
med = np.median(image)
sig = rb.std(image)
vmin = med - siglo*sig
vmax = med + sighi*sig

def loadImage(filename,plot=False):
	image,header = fits.getdata(filename,0,header=True)
	return image,header

def plotImage(filename):
	image = loadImage(filename)
	plt.ion()
	plt.figure(1)
	plt.clf()
	plt.imshow(image,vmin=vmin,vmax=vmax,cmap='gray')

def findSources(filemane,plot=False):
	image,header = loadImage(filename)
	sources = daofind(image - med, fwhm=3.0, threshold=5.*sig)
	positions = (sources['xcentroid'], sources['ycentroid'])
	apertures = CircularAperture(positions, r=4.)
	#norm = ImageNormalize(stretch=SqrtStretch())
	if plot:
		plotImage(filename)
		apertures.plot(color='blue', lw=1.5, alpha=0.5)

	return sources

