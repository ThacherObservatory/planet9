import numpy as np
import matplotlib.pyplot as plt
import planet9_sim as p9s
from astropy.io import fits

oversamp = 10
size = 2048
bias = 500
readnoise = 2
seeing = 3.0
width = 10.0
p9mag = 20.0
exptime = 600.0
p9pos = np.array([1024,1024])*oversamp
plate_scale = 0.61
x, y = p9s.distribute(size=size*oversamp)

image = p9s.make_field(size=size,x=x,y=y,oversamp=oversamp,bias=bias,readnoise=readnoise,seeing=seeing,
		       exptime=exptime,p9pos=p9pos,p9mag=p9mag,plate_scale=plate_scale,width=width)

plt.ion()
plt.figure(1)
plt.clf()
p9s.plot_field(image)
plt.xlim(1000,1048)
plt.ylim(1000,1048)

fits.writeto('20mag_600s.fits', image, clobber=True)
