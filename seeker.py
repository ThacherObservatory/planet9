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
import planet9_sim as p9s
from plot_params import *
import pdb,pickle

#path = "/Users/nickedwards/Desktop/"
#filename = path + "stars.fits"

# some standard variables


def loadImage(filename, plot=True):
	image, header = fits.getdata(filename, 0, header=True)
	med = np.median(image)
	sig = rb.std(image)
	siglo = 3
	sighi = 7
	vmin = med - siglo * sig
	vmax = med + sighi * sig
	if plot:
		plt.ion()
		plt.figure()
		plt.clf()
		plt.imshow(image, vmin=vmin, vmax=vmax,
				   cmap='gray', interpolation='none')

	return image, header, med, sig


def findSources(filename, plot=True, pixscale=0.61, seeing=3.0, threshold=2.0,
				sharplo=0.2, sharphi=1.0, roundlo=-0.5, roundhi=0.5, mzp=22.5):
	image, header, med, sig = loadImage(filename, plot=plot)
	sources = daofind(image - med, fwhm=seeing / pixscale, threshold=threshold * sig,
					  sharplo=sharplo, sharphi=sharphi, roundlo=roundlo, roundhi=roundhi)
	positions = (sources['xcentroid'], sources['ycentroid'])
	apertures = CircularAperture(positions, r=seeing / pixscale)
	if plot:
		apertures.plot(color='cyan', lw=1.5, alpha=0.5)

	return sources


def testRecovery(p9mag=10.0, seeing=3.5, threshold=2.0,
				 sharplo=0.2, sharphi=1.0, roundlo=-0.5, roundhi=0.5,
				 niter=10, debug=False):

	found = 0
	rand = 0
	for i in range(niter):
		p9s.planet9_sequence(readnoise=2, p9mag=p9mag,
							 nimage=1, filename='test')
		sources = findSources('test_1.fits', threshold=threshold, seeing=seeing, sharplo=sharplo,
							  roundlo=roundlo, roundhi=roundhi, plot=False)

		positions = (sources['xcentroid'], sources['ycentroid'])
		apertures = CircularAperture(positions, r=seeing / 0.61)

		d = np.sqrt((1024.0 - sources['xcentroid'])
					** 2 + (1024.0 - sources['ycentroid'])**2)

		xrand, yrand = np.random.uniform(0, 2048, 2)

		drand = np.sqrt(
			(xrand - sources['xcentroid'])**2 + (yrand - sources['ycentroid'])**2)

		# maximum distance is width at half max
		dmax = seeing / (2 * 0.61)

		if np.min(d) < dmax:
                        print 
			found += 1
		if np.min(drand) < dmax:
			rand += 1

		if debug:
			image, header, med, sig = loadImage('test_1.fits', plot=False)

			plt.ion()
			plt.figure(987)
			plt.clf()
		        plt.imshow(image, vmin=med - 2 * sig, vmax=med + 5 *
					   sig, cmap='gray', interpolation='none')
			apertures.plot(color='cyan', lw=1.5, alpha=0.5)
			plt.plot([1024], [1024], 'r+', ms=12)
			plt.xlim(924, 1124)
			plt.ylim(924, 1124)
			plt.show()
			pdb.set_trace()

	percent = np.float(found) / np.float(niter) * 100.0
	control = np.float(rand) / np.float(niter) * 100.0

	return percent, control


def runTest(p9mag=[10, 23], seeing=3.5, threshold=2.0,
			sharplo=0.2, sharphi=1.0, roundlo=-0.5, roundhi=0.5,
			niter=100, debug=False):

	mags = (np.arange(p9mag[0], p9mag[1] + 1)).astype('float')

	percent = []
	control = []
	for mag in mags:
		p, c = testRecovery(p9mag=mag, niter=niter, debug=debug)
		percent = np.append(percent, p)
		control = np.append(control, c)

	plot_params()
	plt.figure(99)
	plt.plot(mags, percent, 'b-')
	plt.plot(mags, control, 'r--')
	plt.xlabel('V Band Magnitude', fontsize=17)
	plt.ylabel('Recovery Probability', fontsize=17)
	plt.savefig('Recovery.png', dpi=300)

        dict = {'mags':mags,'percent':percent,'control':control}
        pickle.dump(dict, open("Recovery.p", "wb"))
        
	return mags, percent, control
