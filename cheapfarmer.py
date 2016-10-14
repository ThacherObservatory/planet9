import matplotlib.pyplot as plt
import numpy as np
import scipy as sp
from tqdm import tqdm
import trylegal as tr
from astropy.io import fits

# Variables
sigma = 20				# sigma
floor = 500				# floor
size = 20480			# image resolution
seeing = 3				# seeing quality
magn_0flux = 22.5		# observatory V-band zero-flux magnitude
plate_scale = 0.61		# plate-scale, arcseconds per pixel
star_width = 10*sigma	# range of star gaussian
exposure_time = 1800	# exposure time in seconds

def make_image():
	# Make an image of given size with specified noise properties.
	# Floor is bias, sigma is noise
	image = np.random.normal(floor,sigma,(size/10,size/10))
	starframe = np.zeros((size,size))

	# Create arrays that correspond to image pixels
	global x
	global y
	x = np.arange(0, size, 1, int)
	y = x[:,np.newaxis]

	return image, starframe

def make_curve():
	x_micro = np.arange(0, (seeing/plate_scale), 1, int)
	y_micro = x[:,np.newaxis]

	x_center = size // 2
	y_center = size // 2

	curve = np.exp(-4*np.log(2) * ((x_micro-x_center)**2 + (y_micro-y_center)**2) / seeing**2)
	return curve

def coord_gen():
	reload(tr)
	rand = []
	pbar = tqdm(desc = 'Randomizing coordinates', total = tr.info_len(), unit = 'whatever(s)')
	for i in range(tr.info_len()):
		value = np.random.randint(0,size)
		rand = np.array(np.append(rand, value))
		pbar.update(1)
	return rand

# generate coordinate sets
def distribute():
	x_rand = coord_gen()
	y_rand = coord_gen()
	return x_rand, y_rand

def star_calc(image, starframe, curve):
	sp.convolve(starframe, curve) # The power of Christ convolves you!
	return image

# place stars into the field
def add_star(starframe, loc, magn):
	# Decide where the star will go
	x0 = int(loc[0][0])
	y0 = int(loc[1][0])
	flux = exposure_time*(10**(-0.4*(magn - magn_0flux))) # flux

	# Add the star into the image
	starframe[x0,y0] = flux
	return starframe

def plot_field(image):
	# Make a plot
	plt.ion()
	plt.figure(1)
	plt.clf()
	plt.imshow(image,cmap='gray')
	plt.gca().invert_yaxis()
	plt.savefig("stars.png")

def slice_plot(image):
	# Show a cross section of the star in the image
	slice = image[size//2,:]
	plt.figure(2)
	plt.clf()
	plt.xlim(size)
	plt.gca().invert_xaxis()
	plt.plot(slice)

def make_field():
	# create coordinate system
	image, starframe = make_image()

	# create random coordinates
	x_rand, y_rand = distribute()
	curve = make_curve()

	tri_data = tr.info_col('V')

	# progress bar
	pbar = tqdm(desc = 'Placing stars', total = tr.info_len(), unit = 'star(s)')

	# generate stars
	for i in range(tr.info_len()):
		loc = [x_rand[[i]], y_rand[[i]]]
		magn = tri_data.iloc[i]['V']
		starframe = add_star(starframe, loc, magn)
		pbar.update(1)
	pbar.close()

	image = star_calc(image, starframe, curve)

	# render image and save it
	fits.writeto('stars.fits', image, clobber = True)
	plot_field(image)
	slice_plot(image)

# run
if __name__ == '__main__':
	make_field()
