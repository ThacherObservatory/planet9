import matplotlib.pyplot as plt
import numpy as np
import scipy as sp
from tqdm import tqdm
import trylegal as tr
from astropy.io import fits
from rebin import rebin
import pdb
import os
import robust as rb

# Variables
# sigma = 20				# Read noise
# floor = 500				# Bias level
# size = 2048				# Image size (pixels)
# oversamp = 10				# Oversampling for subpixel accuracy
# seeing = 3				# Seeing FWHM
# mbg = 21.3				# Background magnitudes/sq. arcsec
# mzp = 22.5	        	# Observatory V-band zeropoint magnitude
# plate_scale = 0.61		# Plate-scale, arcseconds per pixel
# star_width = 10*sigma		# Width of star subframes
# exptime = 1800	        # Exposure time in seconds


def make_noise_frame(size=2048, bias=500, readnoise=20, background=21.3, mzp=22.5, exptime=1800.0,
					 plate_scale=0.61):
	"""
	Make an image of given size with specified noise properties.
	Need to put background noise in here
	"""
	# make image with bias and readnoise
	readimage = np.random.normal(bias, readnoise, (size, size))

	# Get background flux from mags per square arcsec
	bgflux = exptime * (10**(-0.4 * (background - mzp)))
	pixpersqarcsec = plate_scale**2
	bgflux *= pixpersqarcsec

	# make background image
	bgimage = np.random.poisson(bgflux, (size, size))

	# final image is sum of two images
	image = readimage + bgimage

	return image


def make_blank_frame(size=2048, oversamp=1):
	"""
	Make blank frame of proper dimensions
	"""
	return np.zeros((size * oversamp, size * oversamp))


def make_star(seeing=3.0, plate_scale=0.61, width=5.0):
	"""
	Make a gaussian "star" in the center of a frame
	Normalize to a flux of 1
	"""
	swidth_pix = seeing / plate_scale
	size_pix = np.round(swidth_pix * width / 2.0).astype('int')
	x = np.arange(-size_pix, size_pix, 1, int)
	y = x[:, np.newaxis]
	star = np.exp(-4 * np.log(2) * (x**2 + y**2) / swidth_pix**2)
	star /= np.sum(star)

	return star


def coord_gen(size=2048):
	"""
	Generates uniformly random integers between 0 and size
	"""
	reload(tr)
	rand = []
# pbar = tqdm(desc = 'Randomizing coordinates', total = tr.info_len(),
# unit = 'whatever(s)')
	for i in range(tr.info_len()):
		value = np.random.randint(0, size)
	rand = np.array(np.append(rand, value))
#	pbar.update(1)
	return rand


def distribute(size=2048):
	"""
	Create random x and y coordinates (integer values)
	"""
	x_rand = coord_gen(size=size)
	y_rand = coord_gen(size=size)
	return x_rand.astype('int'), y_rand.astype('int')


# place stars into the field
def add_star(starframe, star, loc=[0, 0], mag=0, mzp=22.5, exptime=1800.0):
	# extract x and y values of the star
	x0 = loc[0]
	y0 = loc[1]

	# compute total flux of star
	flux = exptime * (10**(-0.4 * (mag - mzp)))  # flux

	# turn into integers
	starnorm = star * flux
	star_int = np.round(starnorm)

# Do this later
#    # include signal noise
	star_shape = np.shape(star_int)
#    poisson_star = np.random.poisson(star_int.flatten())
#    poisson_star = np.reshape(poisson_star,star_shape)

	# size of star frame
	xs = star_shape[0]
	ys = star_shape[0]

	# size of frame
	xf = np.shape(starframe)[0]
	yf = np.shape(starframe)[1]

	# figure out beginning and ending indices
	startx = max(0, x0 - xs / 2)
	if startx == 0:
		xb = xs / 2 - x0
	else:
		xb = 0
	stopx = min(x0 + xs / 2, xf)
	if stopx == xf:
		xe = stopx - startx
	else:
		xe = xs

	starty = max(0, y0 - ys / 2)
	if starty == 0:
		yb = ys / 2 - y0
	else:
		yb = 0
	stopy = min(y0 + ys / 2, yf)
	if stopy == yf:
		ye = stopy - starty
	else:
		ye = ys

	# Add the star into the image
	try:
		starframe[startx:stopx, starty:stopy] += star_int[xb:xe, yb:ye]
	except:
		print startx, stopx, stopx - startx
		print starty, stopy, stopy - starty
		print xb, xe, xe - xb
		print yb, ye, ye - yb
		pdb.set_trace()
		pass
	return starframe


def plot_field(image, siglo=2.0, sighi=5.0, write=False):
	# Make a plot
		med = np.median(image)
		sig = rb.std(image)
	plt.figure(1)
	plt.clf()
	plt.imshow(image,vmin=med-siglo*sig,vmax=med+sighi*sig,cmap='gray',interpolation='none')
	plt.gca().invert_yaxis()
	if write:
			plt.savefig("stars.png")


def slice_plot(image):
	# Show a cross section of the star in the image
		xsize = np.shape(image)[0]
	slice = image[xsize//2,:]
	plt.figure(2)
	plt.clf()
	plt.xlim(xsize)
	plt.plot(slice)


def make_field(size=2048,x=None,y=None,oversamp=10,bias=500,readnoise=20,seeing=3.0,plate_scale=0.61,width=10.0,
			   background=21.3,mzp=22.5,exptime=1800.0,write=False,
			   p9pos=[1000,1000],p9mag=23.0,plot=False):
	"""
	Make a field of stars with realistic noise properties
	"""

	# create coordinate system
	starframe = make_blank_frame(size=size,oversamp=oversamp)
	noiseframe = make_noise_frame(size=size,bias=bias,readnoise=readnoise,background=background,
								  plate_scale=plate_scale,mzp=mzp,exptime=exptime)
	xs = np.shape(noiseframe)[0]
	ys = np.shape(noiseframe)[1]

	# create random coordinates
	if x == None or y == None:
		x, y = distribute(size=size*oversamp)

	star = make_star(seeing=seeing,plate_scale=plate_scale/oversamp,width=width)

	tri_data = tr.info_col('V')

	# progress bar
	pbar = tqdm(desc = 'Placing stars', total = tr.info_len(), unit = 'stars')

	# generate stars
	for i in range(tr.info_len()):
		loc = [x[i], y[i]]
		mag = tri_data.iloc[i]['V']
	starframe = add_star(starframe, star, loc=loc, mag=mag, mzp=mzp, exptime=exptime)
	pbar.update(1)
	pbar.close()

	# add p9 in
	starframe = add_star(starframe, star, loc=p9pos, mag=p9mag, mzp=mzp, exptime=exptime)

	# rebin oversampled data
	if oversamp > 1:
		starframe = rebin(starframe,xs,ys)

	# add poisson noise to the star image
	shape = np.shape(starframe)
	star_int = starframe.astype('int')
	star_noise = np.random.poisson(star_int.flatten())
	star_noise = np.reshape(star_noise,shape)

	image = star_noise + noiseframe

	# render image and save it
	if write:
		fits.writeto('stars.fits', image, clobber = True)

	if plot:
		plot_field(image,write=write)

	return image


def planet9_movie(size=2048,oversamp=10,bias=500,readnoise=20,seeing=3.0,
					 plate_scale=0.61,width=10.0,background=21.3,mzp=22.5,exptime=1800.0,
					 write=False,p9pos=[1024,1024],p9mag=23.0,dpos=30.0,angle=225.0,nimage=4,
					 filename='P9',fps=2):

	# get locations of stars
	x, y = distribute(size=size*oversamp)

	# get locations of P9
	t = np.arange(0,dpos*nimage,dpos)
	p9_x = p9pos[0]*oversamp + t*np.cos(np.radians(angle))*oversamp/plate_scale
	p9_y = p9pos[1]*oversamp + t*np.sin(np.radians(angle))*oversamp/plate_scale

	for i in range(nimage):
		image = make_field(size=size,x=x,y=y,oversamp=oversamp,bias=bias,readnoise=readnoise,seeing=seeing,
						   plate_scale=plate_scale,width=width,background=background,mzp=mzp,exptime=exptime,
						   write=write,p9pos=[p9_x[i],p9_y[i]],p9mag=p9mag)
		fname = 'p9_image%05d.png'%i
		plt.savefig(fname,bbox_inches='tight',transparent=True, pad_inches=0,frameon=False,
					dpi=150)
		os.system("convert "+fname+" -background black -flatten +matte "+fname)

	os.system("rm "+filename+".mp4")
	os.system("ffmpeg -r "+str(fps)+" -i p9_image%05d.png -b:v 20M -vcodec libx264 -pix_fmt yuv420p -s 808x764 "+\
			  filename+".mp4")
	os.system("rm p9_image*png")


def planet9_sequence(size=2048,oversamp=10,bias=500,readnoise=20,seeing=3.0,
					 plate_scale=0.61,width=10.0,background=21.3,mzp=22.5,exptime=1800.0,
					 write=False,p9pos=[1024,1024],p9mag=23.0,dpos=30.0,angle=225.0,nimage=4,
					 filename='P9'):
	# get locations of stars
	x, y = distribute(size=size*oversamp)

	# get locations of P9
	t = np.arange(0,dpos*nimage,dpos)
	p9_x = p9pos[0]*oversamp + t*np.cos(np.radians(angle))*oversamp/plate_scale
	p9_y = p9pos[1]*oversamp + t*np.sin(np.radians(angle))*oversamp/plate_scale

	for i in range(nimage):
		image = make_field(size=size,x=x,y=y,oversamp=oversamp,bias=bias,readnoise=readnoise,seeing=seeing,
						   plate_scale=plate_scale,width=width,background=background,mzp=mzp,exptime=exptime,
						   write=write,p9pos=[p9_x[i],p9_y[i]],p9mag=p9mag)

		# write image to FITS
		fits.writeto(filename+'_'+str(i+1)+'.fits', image, clobber = True)
	return



# run
if __name__ == '__main__':
	planet9_movie()
