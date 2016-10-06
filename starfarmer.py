import matplotlib.pyplot as plt
import numpy as np
from tqdm import tqdm
import trylegal as tr

# Variables
snr = 20		# signal-to-noise ratio
sigma = 20		# sigma
floor = 500		# floor
size = 2048		# image resolution
seeing = 3		# seeing quality
flux = 100		# flux

def make_image():
	# Make an image of given size with specified noise properties.
	# Floor is bias, sigma is noise
	image = np.random.normal(floor,sigma,(size,size))

	# Create arrays that correspond to image pixels
	x = np.arange(0, size, 1, float)
	y = x[:,np.newaxis]
	return image, x, y

def coord_gen():
	reload(tr)
	rand = []
	for i in range(tr.info_len()):
		value = np.random.uniform(0,size)
		rand = np.array(np.append(rand, value))
	return rand

# generate coordinate sets
def distribute():
	x_rand = coord_gen()
	y_rand = coord_gen()
	return x_rand, y_rand

# place stars into the field
def add_star(image, x, y, loc, snr, sigma, floor, seeing, flux):
	# Interpret flux in terms of noise level
	flux *= sigma

	# Decide where the star will go. If no input, place the star in the
	# middle of the image
	if loc is None:
		x0 = y0 = size // 2
	else:
		x0 = loc[0]
		y0 = loc[1]

	# This turns the star into a Gaussian
	star = np.exp(-4*np.log(2) * ((x-x0)**2 + (y-y0)**2) / seeing**2)

	# This is approximately how many pixels the star covers (needed
	# to compute the total signal to noise)
	npix = max(np.pi*(seeing)**2,1)

	# Total noise within the star image
	noise = np.sqrt(npix)*sigma

	# Adjust brightness of the star if SNR is specified in the inputs
	if snr:
		bright = noise*snr
	elif flux:
		bright = flux
		snr = flux/noise
	star  = star/np.sum(star) * bright

	# Add the star into the image
	image += star

	return image

def plot_field(image):
	# Make a plot
	plt.ion()
	plt.figure(1)
	plt.clf()
	plt.imshow(image,cmap='gray')
	plt.savefig("stars.png")

def slice_plot(image,):
	# Show a cross section of the star in the image
	slice = image[size//2,:]
	plt.figure(2)
	plt.clf()
	plt.xlim(size)
	plt.gca().invert_xaxis()
	plt.plot(slice)

def make_field():
	# create coordinate system
	image, x, y = make_image()

	# create random coordinates
	x_rand, y_rand = distribute()

	# progress bar
	pbar = tqdm(desc = 'Creating image', total = tr.info_len(), leave = False, unit = 'star(s)')

	# generate stars
	for i in range(tr.info_len()):
		loc = [x_rand[[i]], y_rand[[i]]]
		add_star(image, x, y, loc, snr, sigma, floor, seeing, flux)
		pbar.update(1)
	pbar.close()

	# render image
	plot_field(image)
	slice_plot(image)

# run
if __name__ == '__main__':
	make_field()
