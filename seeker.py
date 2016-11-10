import matplotlib.pyplot as plt
import numpy as np
import scipy as sp
from astropy.io import fits
from astropy.stats import sigma_clipped_stats
# talk to doc swift about this
#from photutils import DAOStarFinder
from photutils import daofind
from photutils import CircularAperture
from astropy.visualization import SqrtStretch
from astropy.visualization.mpl_normalize import ImageNormalize
import robust as rb


path = "/Users/nickedwards/Desktop/"

filename = path + "stars.fits"

# some standard variables



def loadImage(filename,plot=True):
	image,header = fits.getdata(filename,0,header=True)
	med = np.median(image)
	sig = rb.std(image)
	siglo=3
	sighi=7
	vmin = med - siglo*sig
	vmax = med + sighi*sig

	if plot:
		plt.ion()
		plt.figure()
		plt.clf()
		plt.imshow(image,vmin=vmin,vmax=vmax,cmap='gray')

	return image,header,med,sig

def findSources(filename,plot=True):
	image,header,med,sig = loadImage(filename,plot=plot)
	sources = daofind(image - med, fwhm=3.0/.61, threshold=2*sig,sharphi=.6)
	positions = (sources['xcentroid'], sources['ycentroid'])
	apertures = CircularAperture(positions, r=4.)
	#norm = ImageNormalize(stretch=SqrtStretch())
	if plot:
		apertures.plot(color='blue', lw=1.5, alpha=0.5)

	return sources

