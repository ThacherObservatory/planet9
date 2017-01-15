import numpy as np
import matplotlib.pyplot as plt
from quick_image import display_image, readimage
import thacherphot as tp
import hcongrid as h
from astropy.io import fits
from kapteyn import maputils
import sys

bias = tp.master_bias(None)
flat = tp.master_flat(None)

datadir = '/Users/jonswift/Astronomy/ThacherObservatory/data/26Dec2016/'
files,fct = tp.get_files(dir=datadir,prefix='Orion',suffix='solved.fits')

zsz = len(files)
reffile = files[zsz/2]

image0,header0 = readimage(reffile)

ysz, xsz = np.shape(image0)
refim = h.pyfits.open(reffile)
refh = h.pyfits.getheader(reffile)
stack = np.zeros((xsz,ysz,zsz))

for i in range(zsz):
    print 'Starting image '+str(i)+' of '+str(zsz)
    im = h.pyfits.open(files[i])
    newim = h.hcongrid((im[0].data-bias)/flat, im[0].header,refh)
    stack[:,:,i] = newim
    
final = np.median(stack, axis=2)

display_image(final)

fits.writeto('P9_sample_image.fits', final, refh)

sys.exit()


from kapteyn import maputils
from matplotlib import pyplot as plt
import numpy as np
from quick_image import readimage

image0,header0 = readimage("P9_sample_image.fits")

clipmin = np.median(image0)-0.3*np.std(image0)
clipmax = np.median(image0)+2*np.std(image0)


f = maputils.FITSimage("P9_sample_image.fits")
fig = plt.figure()
frame = fig.add_subplot(1,1,1)
annim = f.Annotatedimage(frame,cmap="gray", clipmin=clipmin, clipmax=clipmax)
annim.Image()
grat = annim.Graticule()
grat.setp_gratline(visible=True,linestyle='--')
grat.setp_axislabel("bottom",label=r"$\mathrm{Right\ Ascension\ (J2000)}$",fontsize=16)
grat.setp_axislabel("left",label=r"$\mathrm{Declination\ (J2000)}$",fontsize=16)
grat.setp_ticklabel(plotaxis='left',fontsize=14)
grat.setp_ticklabel(plotaxis='bottom',fontsize=14)
annim.plot()
plt.savefig('P9_sample_image.png',dpi=300)
plt.show()

