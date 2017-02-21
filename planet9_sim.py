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


"""
==================
Variables Manifest
==================

* denotes variables that are manually defined
^ denotes optional switches

*angle (float):         planet 9's direction of travel from starting position, degrees
*background (float):    ???
*bias (float):          camera bias level
*dpos (float):          planet 9's change in position between frames
*exptime (float):       exposure time in seconds
*filename (str):        filename for saving video file
*fps (float):           frames per second for the video
*gain (float):          camera gain level
^grid (bool):           whether or not to overlay a grid
image (array):          array to store simulated image data
loc (array):            coordinate pair for the simulated star
mag (float):            apparent magnitude of the simulated star
*mzp (float):           mean zero-point value
*nimage (int):          number of frames to generate
*oversamp (int):        scaling factor for rendering the simulated data
*p9mag (float):         simulated planet 9 apparent magnitude
*p9pos (array):         starting coordinate pair for simulated planet 9
*plate_scale (float):   width of 1 camera pixel in arcseconds
^plot (bool):           whether or not to plot the simulated image
*readnoise (float):     camera read noise level
*seeing (float):        seeing quality in arcseconds
*sighi (float):         ???
*siglo (float):         ???
*size (int):            image height and width in pixels
star (array):           array of simulated star data
starframe (array):      array to store simulated image data
*stretch (str):         stretch method to apply to simulated star field for viewing
*width (float):         ???
^write (bool):          whether or not to save the chart as an image file
^write (bool):          whether or not to save to an image file
^write (bool):          whether or not to write the video to a file
"""

# Image frame settings
size = int(2048)

# Rendering settings
oversamp = 10
stretch = 'linear'
siglo = -1.0
sighi = 5.0

# Virtual camera settings
bias = 500
exptime = 1800.0
gain = 6.595
plate_scale = 0.61
readnoise = 20

# Virtual observational variables
background = 21.3
mzp = 22.5
seeing = 3.0
width = 10.0

# Planet 9 variables
angle = 225.0
dpos = 30.0
p9mag = 22.0
p9pos = [1024, 1024]


def make_noise_frame(background):
    """Make an image of given size with specified noise properties.
    Need to put background noise in here

    Variables:
        readimage (array):
        bgflux (float):
        pixpersqarcsec (float):
        ecounts (float):
        bgimage (array):

    Returns:
        noise_frame (array):          composite array of simulated readings and simulated noise
    """
    # make image with bias and readnoise
    readimage = np.random.normal(bias, readnoise, (size, size))

    # Get background flux from mags per square arcsec
    bgflux = exptime * (10**(-0.4 * (background - mzp)))
    pixpersqarcsec = plate_scale**2
    bgflux *= pixpersqarcsec

    ecounts = bgflux * gain

    # make background image
    bgimage = (np.random.poisson(ecounts, (size, size))) / gain

    # final image is sum of two images
    noise_frame = readimage + bgimage

    return noise_frame


def make_blank_frame(oversamp):
    """Make blank frame of proper dimensions

    Args:
        oversamp (int):         factor by which to scale the image

    Returns:
        blank frame (array):    empty array for the image
    """
    blankframe = np.zeros((size * oversamp, size * oversamp))
    return blankframe


def make_source_frame():
    """Make a gaussian "star" in the center of a frame
    Normalize to a flux of 1

    Args:
        width (float):          width of area to calculate to render star

    Variables:
        swidth_pix (float):
        size_pix (float):
        x (array):              coordinates for the area to calculate the source
        y (array):              "

    Returns:
        star (array):           array with simulated data for a single star
    """
    swidth_pix = seeing / plate_scale
    size_pix = np.round(swidth_pix * width / 2.0).astype('int')
    x = np.arange(-size_pix, size_pix, 1, int)
    y = x[:, np.newaxis]
    star = np.exp(-4 * np.log(2) * (x**2 + y**2) / swidth_pix**2)
    star /= np.sum(star)

    return star


def coord_gen(coord_limit):
    """Generates uniformly random integers between 0 and size

    Arguments:
        coord_limit (int):      upper bourndary for random coordinate generation

    Variables:
        value (int):            randomly generated number for rand[]

    Returns:
        rand (array):           list of random values for each simulated star
    """
    rand = []
    for i in range(tr.info_len()):
        value = np.random.randint(0, coord_limit)
        rand = np.append(rand, value)
    return rand


def distribute(coord_limit):
    """Create random x and y coordinates (integer values)

    Arguments:
        coord_limit (int):      upper bourndary for random coordinate generation

    Returns:
        x_rand (array):         list of random values to be used as x-coordinates
        y_rand (array):         list of random values to be used as y-coordinates
    """
    x_rand = coord_gen(coord_limit)
    y_rand = coord_gen(coord_limit)
    return x_rand, y_rand


def distribute_oversamp():
    """Expands generated star coordinate-space for oversampling

    Variables:
        x (array):
        y (array):
    """

    # get locations of stars
    x, y = distribute(size * oversamp)

    return x, y


def add_source(starframe, sourceframe, loc, mag):
    """places sources into the image frame array

    Args:
        starframe (array):      array to store simulated image data
        sourceframe (array):    array of simulated source data
        loc (array):            coordinate pair for the simulated star
        mag (float):            apparent magnitude of the simulated star

    Variables:
        x0 (int):
        y0 (int):
        flux (float):
        starnorm (array):       generic star curve adjusted for individual flux
        star_int (array):       starnorm, reduced to integers and 1 dimension
        star_shape (tuple):     dimensions of starnorm
        xs (???):
        ys (???):
        xf (???):
        yf (???):
        xb (???):
        yb (???):
        xe (???):
        ye (???):
        startx (int):
        stopx (int):
        starty (int):
        stopy (int):

    Returns:
        starframe (array):      array to store simulated image data, with star[] added
    """
    # extract x and y values of the star
    x0 = loc[0]
    y0 = loc[1]

    # compute total flux of star
    # flux = exptime * (10**(-0.4 * (mag - mzp))) * (oversamp**2)  # flux
    flux = exptime * (10**(-0.4 * (mag - mzp)))

    # turn into integers
    starnorm = sourceframe * flux
    star_int = np.round(starnorm)

    # Do this later
    # include signal noise
    # poisson_star = np.random.poisson(star_int.flatten())
    # poisson_star = np.reshape(poisson_star,np.shape(star_int))

    # size of star frame
    xs = np.shape(star_int)[0]
    ys = np.shape(star_int)[0]

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
        starframe[np.int(startx):np.int(stopx), np.int(starty):np.int(stopy)] \
            += star_int[np.int(xb):np.int(xe), np.int(yb):np.int(ye)]
    except:
        print startx, stopx, stopx - startx
        print starty, stopy, stopy - starty
        print xb, xe, xe - xb
        print yb, ye, ye - yb
        pdb.set_trace()
        pass
    return starframe


def plot_field(image, grid=False, write=False):
    """generates a plot to display the simulated field image

    Args:
        image (array):          array to store simulated image data
        siglo (float):          ???
        sighi (float):          ???
        grid (bool):            whether or not to overlay a grid
        write (bool):           whether or not to save the chart as an image file
    """
    # Make a plot
    if stretch == 'sqrt':
        image = np.sqrt(image)
    med = np.median(image)
    sig = rb.std(image)
    plt.ion()
    plt.figure(1)
    plt.clf()
    plt.imshow(image, vmin=med - siglo * sig, vmax=med +
               sighi * sig, cmap='gray', interpolation='none')
    if grid:
        plt.rc('grid', linestyle='-', color='white')
        plt.grid(which='both')
    plt.gca().invert_yaxis()
    if write:
        plt.savefig("stars.png", dpi=300)


def slice_plot(image):
    """Show a cross section of the star in the image

    Args:
        image (array):          array to store simulated image data
    """
    xsize = np.shape(image)[0]
    slice = image[xsize // 2, :]
    plt.ion()
    plt.figure(2)
    plt.clf()
    plt.xlim(xsize)
    plt.plot(slice)


def make_field(x, y, background, p9pos, plot=False, grid=False, write=False):
    """Make a field of stars with realistic noise properties

    Args:
        p9pos (array):          starting coordinate pair for simulated planet 9
        plot (bool):            whether or not to plot the simulated image
        grid (bool):            whether or not to overlay a grid
        write (bool):           whether or not to save to an image file

    Variables:
        noiseframe (array):
        star (array):
        loc (array):
        mag (float):
        starframe (array):
        shape (tuple):
        star_int (tuple):
        star_noise (???):
        star_noise (???):
        tri_data (dataframe):
        xs (array):
        ys (array):

    Returns:
        image (array):          array containing the simulated image data
    """

    # create coordinate system
    starframe = make_blank_frame(oversamp)
    noiseframe = make_noise_frame(background)
    xs = np.shape(noiseframe)[0]
    ys = np.shape(noiseframe)[1]

    # do something if no coordinates

    sourceframe = make_source_frame()

    tri_data = tr.info_col('V')

    # progress bar
    pbar = tqdm(desc='Placing stars', total=tr.info_len(), unit='stars')

    # generate stars
    for i in range(tr.info_len()):
        loc = [x[i], y[i]]
        mag = tri_data.iloc[i]['V']
        starframe = add_source(starframe, sourceframe, loc, mag)
        pbar.update(1)
    pbar.close()

    # add p9 in
    starframe = add_source(starframe, sourceframe, p9pos, p9mag)

    # rebin oversampled data
    if oversamp > 1:
        starframe = rebin(starframe, xs, ys)

    # add poisson noise to the star image
    shape = np.shape(starframe)
    stars_int = np.round(starframe).astype('int')
    stars_noise = np.random.poisson(stars_int.flatten())
    stars_noise = np.reshape(stars_noise, shape)

    # print "Verify poission noise, line 393"
    # pdb.set_trace()
    image = stars_noise + noiseframe
    # render image and save it
    if write:
        fits.writeto('stars.fits', image, clobber=True)

    if plot:
        plot_field(image, grid, write)

    return image
    # return np.sqrt(image)


def planet9_path(nimage):
    """Generate planet 9 coordinates

    Variables:
        t (array):
        p9_x (array):
        p9_y (array):
    """
    # get locations of P9, move it for each frame
    t = np.arange(0, dpos * nimage, dpos)
    p9_x = p9pos[0] * oversamp + t * np.cos(np.radians(angle)) * oversamp / plate_scale
    p9_y = p9pos[1] * oversamp + t * np.sin(np.radians(angle)) * oversamp / plate_scale

    return p9_x, p9_y


def planet9_movie(nimage=4, fps=2, grid=False, write=False, filename='P9'):
    """create a video from simulated field images

    Args:
        nimage (int):           number of frames to generate
        fps (float):            frames per second for the video
        grid (bool):            whether or not to overlay a grid
        write (bool):           whether or not to write the video to a file
        filename (str):         filename for saving video file

    Variables:
        x (array):
        y (array):
        t (array):
        p9_x (array):
        p9_y (array):
        p9pos (array):          current coordinate pair for planet 9
        image (array):
        fname (str):
    """

    # initialize persistent data
    p9_x, p9_y = planet9_path(nimage)
    x, y = distribute_oversamp()

    # proress bar
    pbar = tqdm(desc='Rendering frames', total=nimage, unit='frame')

    for i in range(nimage):
        p9pos = [p9_x[i], p9_y[i]]
        image = make_field(x, y, background, p9pos)

        fname = 'p9_image%05d.png' % i

        # grid attempt
        plot_field(image, grid, write)
        plt.savefig(fname, bbox_inches='tight', transparent=True, pad_inches=0, frameon=False,
                    dpi=150)
        os.system("convert " + fname + " -background black -flatten +matte " + fname)
        pbar.update(1)
    pbar.close()
    os.system("rm " + filename + ".mp4")
    os.system("ffmpeg -r " + str(fps) + " -i p9_image%05d.png -b:v 20M -vcodec libx264 -pix_fmt yuv420p -s 808x764 " +
              filename + ".mp4")
    os.system("rm p9_image*png")


# run
if __name__ == '__main__':
    planet9_movie()
