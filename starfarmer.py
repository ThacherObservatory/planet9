import matplotlib.pyplot as plt
import numpy as np
from tqdm import tqdm
import trylegal as tr
from astropy.io import fits

# Variables
sigma = 20				# sigma
floor = 500				# floor
size = 2048				# image resolution
seeing = 3				# seeing quality
magn_0flux = 22.5		# observatory V-band zero-flux magnitude
plate_scale = 0.61		# plate-scale, arcseconds per pixel
star_width = 10*sigma	# range of star gaussian
exposure_time = 1800	# exposure time in seconds

def make_image():
	# Make an image of given size with specified noise properties.
	# Floor is bias, sigma is noise
	image = np.random.normal(floor,sigma,(size,size))

	# Create arrays that correspond to image pixels
	global x
	global y
	x = np.arange(0, size, 1, float)
	y = x[:,np.newaxis]

	return image

def coord_gen():
	reload(tr)
	rand = []
	pbar = tqdm(desc = 'Randomizing coordinates', total = tr.info_len(), unit = 'whatever(s)')
	for i in range(tr.info_len()):
		value = np.random.uniform(0,size)
		rand = np.array(np.append(rand, value))
		pbar.update(1)
	return rand

# generate coordinate sets
def distribute():
	x_rand = coord_gen()
	y_rand = coord_gen()
	return x_rand, y_rand

def star_calc(loc, magn):
	# Interpret flux in terms of noise level
	flux = exposure_time*(10**(-0.4*(magn - magn_0flux))) # flux
	flux *= sigma

	# Decide where the star will go
	x0 = loc[0]
	y0 = loc[1]

	# This turns the star into a Gaussian, only in a star_width*star_width square
	#star = np.exp(-4*np.log(2) * (((x[x0-(star_width/2):x0+(star_width/2),])-x0)**2 + ((y[y0-(star_width/2):y0+(star_width/2),])-y0)**2) / seeing**2)
	star = np.exp(-4*np.log(2) * ((x-x0)**2 + (y-y0)**2) / seeing**2)

	star = star/np.sum(star) * flux

	# Poisson noise
	star = np.random.poisson(star)

	return star

# place stars into the field
def add_star(image, loc, magn):
	star = star_calc(loc, magn)

	# Decide where the star will go
	x0 = loc[0]
	y0 = loc[1]

	# Add the star into the image
	target = image[x0-(star.shape[0]/2):x0+(star.shape[0]/2), y0-(star.shape[1]/2):y0+(star.shape[1]/2)]

	image += star[(target.shape[0]), (target.shape[1]),]
	return image

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
	image = make_image()

	# create random coordinates
	x_rand, y_rand = distribute()
	tri_data = tr.info_col('V')

	# progress bar
	pbar = tqdm(desc = 'Creating image', total = tr.info_len(), unit = 'star(s)')

	# generate stars
	for i in range(tr.info_len()):
		loc = [x_rand[[i]], y_rand[[i]]]
		magn = tri_data.iloc[i]['V']
		image = add_star(image, loc, magn)
		pbar.update(1)
	pbar.close()

	# render image and save it
	fits.writeto('stars.fits', image, clobber = True)
	plot_field(image)
	slice_plot(image)

# run
if __name__ == '__main__':
	make_field()
