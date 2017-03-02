#!/usr/bin/env python
import matplotlib.pyplot as plt
import numpy as np
import scipy as sp
from tqdm import tqdm
import trylegal as tr
from astropy.io import fits
from rebin import rebin
import pdb
import os
import sys
import subprocess
from subprocess import call
import gc
import glob
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
width = 5

# Planet 9 variables
angle = 225.0
dpos = 30.0
p9mag = 22.0
p9pos = [1024, 1024]


def load_data():
    """Initialize data and variables for the script.
    Loads TRILEGAL data

    Returns:
        data (df):              Pandas dataframe of star data
    """
    # load TRILEGAL data
    data = tr.info_col('V')

    return data


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
        source_res (float):       width in pixels of the simulated star
        frame_res (float):      width in pixels of the field to calculate the gaussian
        x (array):              coordinates for the area to calculate the source
        y (array):              "

    Returns:
        star (array):           array with simulated data for a single star
    """
    source_res = (seeing / plate_scale) * oversamp
    frame_res = np.round(source_res * width / 2.0).astype('int')
    x = np.arange(-frame_res, frame_res, 1, int)
    y = x[:, np.newaxis]
    sourceframe = np.exp(-4 * np.log(2) * (x**2 + y**2) / source_res**2)
    sourceframe /= np.sum(sourceframe)

    return sourceframe


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
        x (array):              x-coordinates for simulated stars
        y (array):              y-coordinates for simulated stars
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
        source_int (array):       starnorm, reduced to integers and 1 dimension
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
    flux = exptime * (10**(-0.4 * (mag - mzp))) * (oversamp**2)  # flux

    # turn into integers
    source_convolved = sourceframe * flux
    source_convolved = rebin(source_convolved, np.shape(
        source_convolved)[0], np.shape(source_convolved)[1])
    source_int = np.round(source_convolved).astype(int)

    # Do this later
    # include signal noise
    # poisson_star = np.random.poisson(source_int.flatten())
    # poisson_star = np.reshape(poisson_star,np.shape(source_int))

    # size of star frame
    xs = np.shape(source_int)[0]
    ys = np.shape(source_int)[0]

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
        starframe[np.round(startx).astype(int):np.round(stopx).astype(int), np.round(starty).astype(int):np.round(stopy).astype(int)] \
            += source_int[np.round(xb).astype(int):np.round(xe).astype(int), np.round(yb).astype(int):np.round(ye).astype(int)]
    except:
        print(startx, stopx, stopx - startx)
        print(starty, stopy, stopy - starty)
        print(xb, xe, xe - xb)
        print(yb, ye, ye - yb)
        pdb.set_trace()
        pass
    return starframe


def make_field(source_data, x, y, background, p9pos, oversamp, write, plot=False, grid=False):
    """Make a field of stars with realistic noise properties

    Args:
        source_data
        x
        y
        background
        p9pos (array):          starting coordinate pair for simulated planet 9
        oversamp
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
        source_int (tuple):
        star_noise (???):
        star_noise (???):
        tri_data (dataframe):
        pbar (tqdm):            progress bar for generating the frame
        sim_pbar (tqdm):        progress bar for simulating data

    Returns:
        image (array):          array containing the simulated image data
    """

    # progress bar for the frame
    tasks = 9
    if plot:
        tasks += 1
    if grid:
        tasks += 1
    if write:
        tasks += 1
    pbar = tqdm(desc='Rendering frame', total=tasks, unit='operation', leave=False)

    # create frame
    starframe = make_blank_frame(oversamp)
    pbar.update(1)

    # create convolution frame
    sourceframe = make_source_frame()
    pbar.update(1)
    pbar.refresh()

    # generate stars
    pbar.write('		Simulating stars...')
    sim_pbar = tqdm(desc='Simulating stars', total=tr.info_len(), unit='star', leave=False)
    for i in range(tr.info_len()):
        loc = [x[i], y[i]]
        mag = source_data.iloc[i]['V']
        starframe = add_source(starframe, sourceframe, loc, mag)
        sim_pbar.update(1)
    sim_pbar.close()
    pbar.write('		Stars simulated.')
    pbar.update(1)

    # add p9 in
    starframe = add_source(starframe, sourceframe, p9pos, p9mag)
    pbar.update(1)

    # create noise frame
    noiseframe = make_noise_frame(background)
    pbar.update(1)

    # rebin oversampled data
    if oversamp > 1:
        starframe = rebin(starframe, np.shape(noiseframe)[0], np.shape(noiseframe)[1])
        pbar.update(1)

    # convert data to integers
    stars_int = np.round(starframe).astype('int')
    pbar.update(1)

    # add poisson noise to the star image
    stars_noise = np.random.poisson(stars_int.flatten())
    stars_noise = np.reshape(stars_noise, np.shape(starframe))
    pbar.update(1)

    # print("Verify poission noise, line 393")
    # pdb.set_trace()
    image = stars_noise + noiseframe
    pbar.update(1)

    # render image and save it
    if write:
        fname = 'p9_frame%02d.fits' % it_num
        fits.writeto(fname, image, overwrite=True)
        pbar.update(1)
        pbar.write('		Frame saved as FITS.')

    if plot:
        plot_field(image, grid, write)
        pbar.update(1)

    pbar.close()

    return image


def plot_field(image, grid=False, write=False):
    """generates a plot to display the simulated field image

    Args:
        image (array):          array to store simulated image data
        siglo (float):          ???
        sighi (float):          ???
        grid (bool):            whether or not to overlay a grid
        write (bool):           whether or not to save the chart as an image file
    """
    # progress bar
    tasks = 2
    if grid:
        tasks += 1
    if write:
        tasks += 1
    pbar = tqdm(desc='Plotting image', total=tasks, unit='operation', leave=False)
    pbar.write('		Plotting frame...')

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
    pbar.update(1)

    if grid:
        plt.rc('grid', linestyle='-', color='white')
        plt.grid(which='both')
        pbar.update(1)
    plt.gca().invert_yaxis()
    pbar.update(1)
    if write:
        fname = 'p9_stars%02d.png' % it_num
        plt.savefig(fname, dpi=300)
        pbar.update(1)
        pbar.write('		Frame saved as image.')
    pbar.write('		Plot rendered.')
    pbar.close()


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


def planet9_path(nimage, oversamp):
    """Generate planet 9 coordinates

    Args:
        nimage
        oversamp

    Variables:
        t (array):              array of position for planet 9 throughout frames
        p9_x (array):           set of x-coordinates for planet 9
        p9_y (array):           set of y-coordinates for planet 9
    """
    # get locations of P9, move it for each frame
    t = np.arange(0, dpos * nimage, dpos)
    p9_x = p9pos[0] * oversamp + t * np.cos(np.radians(angle)) * oversamp / plate_scale
    p9_y = p9pos[1] * oversamp + t * np.sin(np.radians(angle)) * oversamp / plate_scale

    return p9_x, p9_y


def frame_render(source_data, x, y, p9pos, oversamp, grid, write, pbar):
    """Render a single frame of simulated noise and data, plus planet 9 and save it

    Args:
        source_data
        x
        y
        p9pos
        i
        pbar (tqdm):            progress bar

    Variables:
        image (array):          array of the simulated data
        fname (str):            filename format for image frames
    """
    pbar.write('	Rendering frame %d...' % (it_num + 1))
    image = make_field(source_data, x, y, background, p9pos, oversamp, write)

    fname = 'p9_frame%02d.png' % it_num

    # grid attempt
    plot_field(image, grid, write)
    plt.savefig(fname, bbox_inches='tight', transparent=True, pad_inches=0, frameon=False,
                dpi=150)
    convert = subprocess.Popen(["magick", "convert", fname, "-background", "black",
                                "-flatten", "+matte", fname], stdout=subprocess.PIPE)
    if convert.communicate()[0]:
        pbar.write(convert.communicate()[0])
    pbar.update(1)
    pbar.write('	Frame %d rendered.' % (it_num + 1))


def planet9_movie(nimage=4, fps=2, grid=False, write=True, filename='P9'):
    """create a video from simulated field images

    Args:
        nimage (int):           number of frames to generate
        fps (float):            frames per second for the video
        grid (bool):            whether or not to overlay a grid
        write (bool):           whether or not to write the video to a file
        filename (str):         filename for saving video file

    Variables:
        x (array):              oversampled x-coordinates for simulated image
        y (array):              oversampled y-coordinates for simulated image
        p9_x (array):           set of x-coordinates for planet 9
        p9_y (array):           set of y-coordinates for planet 9
        p9pos (array):          current coordinate pair for planet 9
        pbar (tqdm):            progress bar for overall process
        frames_pbar (tqdm):     progress bar for generating frames
    """
    # proress bar
    tasks = nimage + 4
    if os.path.exists(filename + ".mp4"):
        tasks += 1
    pbar = tqdm(desc='Simulation progress', total=tasks, unit='operation')

    # initialize persistent data
    pbar.write('Loading TRILEGAL data...')
    source_data = load_data()
    pbar.write('TRILEGAL data loaded.')
    pbar.update(1)

    p9_x, p9_y = planet9_path(nimage, oversamp)
    x, y = distribute_oversamp()
    pbar.update(1)

    # render frames
    frames_pbar = tqdm(desc='Iterating frames', total=nimage, unit='operation', leave=False)
    for i in range(nimage):
        global it_num
        it_num = i
        frame_render(source_data, x, y, [p9_x[i], p9_y[i]], oversamp, grid, write, frames_pbar)
        pbar.update(1)
    frames_pbar.close()

    # remove previous movie
    if os.path.exists(filename + ".mp4"):
        os.remove(filename + ".mp4")
        pbar.update(1)
        pbar.write('Previous movie file deleted.')

    # export new movie with FFmpeg
    ffmpeg = subprocess.Popen(["ffmpeg", "-loglevel", "quiet", "-r", str(fps), "-i", "p9_frame%02d.png",
                               "-b:v", "20M", "-vcodec", "libx264", "-pix_fmt", "yuv420p", "-s", "808x764", filename + ".mp4"], stdout=subprocess.PIPE)

    if ffmpeg.communicate()[0]:
        pbar.write(ffmpeg.communicate()[0])
    else:
        pbar.write('Movie exported.')
        pbar.update(1)

    # delete frame images
    for f in glob.glob("p9_frame??.png"):
        os.remove(f)
    pbar.write('Frame images deleted.')
    pbar.update(1)

    pbar.close()
    gc.collect()
    sys.modules[__name__].__dict__.clear()


# run
if __name__ == '__main__':
    print('Beginning simulation...')
    planet9_movie()
    print('Simulation complete.')
