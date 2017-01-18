#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun Jan 15 13:44:34 2017

@author: ONeill
"""

import numpy as np
import matplotlib.pyplot as plt
from Quickimage import show_image, readimage
import thacherphot as tp
import hcongrid as h
from astropy.io import fits
from kapteyn import maputils
import sys

datadir = '/Users/ONeill/astronomy/data/14Jan2017/'

bias_files = tp.get_files(dir=datadir,prefix='P9',suffix='bias.fit')
flat_files = tp.get_files(dir=datadir,prefix='cal',suffix='.fit')
bias = tp.master_bias(bias_files)
flat = tp.master_flat(flat_files)

files,fct = tp.get_files(dir=datadir,prefix='P9',suffix='solved.fits')

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

show_image(final)

fits.writeto('P9_sample_image.fits', final, refh)

sys.exit()

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

