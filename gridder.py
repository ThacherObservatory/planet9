# -*- coding: utf-8 -*-
"""
Created on Tue Jan 31 09:46:15 2017

@author: Liam
"""

import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord
import pandas as pd

def grid(num):

    #coordinate for center of field
    center = SkyCoord('5 35 17.3	-5 23 28', unit=(u.hourangle, u.deg))

    #camera frame width (arcminutes)
    camera_frame = 21

    #desired overlap (arcminutes)
    overlap = 1

    shift = camera_frame - overlap
    center_ra = (center.ra.hour)
    center_dec = center.dec.degree

    dy = (shift/60.0)*u.degree # change in dec
    dx = (shift/900.0)*u.hour # change in ra

    delta = spiral(num,num)
    print 'delta'

    coords = np.array([])
    print 'coords'
    for i in range(num**2):
        coords = np.append(coords, SkyCoord(ra = center_ra*u.hour+delta[i,0]*dx, dec = center_dec*u.degree+delta[i,1]*dy, frame = 'icrs'))

    print 'SkyCoord'
    RA = [coords[i].ra.hour for i in range(len(coords))]
    DEC = [coords[i].dec.degree for i in range(len(coords))]
    RA = np.array(RA)
    DEC = np.array(DEC)
    coords = pd.DataFrame({'ra': RA, 'dec': DEC})

    return coords
def gridOld():

    #coordinate for center of field
    center = SkyCoord('5 35 17.3	-5 23 28', unit=(u.hourangle, u.deg))

    #camera frame width (arcminutes)
    camera_frame = 21

    #desired overlap (arcminutes)
    overlap = 1

    shift = camera_frame - overlap
    center_ra = (center.ra.hour)
    center_dec = center.dec.degree

    up1 = SkyCoord(ra=center_ra*u.hour, dec=(center_dec+((shift)/60.0))*u.degree, frame='icrs')
    up1_left1 = SkyCoord(ra=(center_ra+((shift)*1.0/900.0/(np.cos(np.degrees(up1.dec)))))*u.hour, dec=up1.dec, frame='icrs')
    up1_right1 = SkyCoord(ra=(center_ra-((shift)*1.0/900.0/(np.cos(np.degrees(up1.dec)))))*u.hour, dec=up1.dec, frame='icrs')

    left1 = SkyCoord(ra=(center_ra+((shift)*1.0/900.0/(np.cos(np.degrees(center_dec)))))*u.hour, dec=up1.dec, frame='icrs')
    right1 = SkyCoord(ra=(center_ra-((shift)*1.0/900.0/(np.cos(np.degrees(center_dec)))))*u.hour, dec=up1.dec, frame='icrs')

    down1 = SkyCoord(ra=center_ra*u.hour, dec=(center_dec-((shift)/60.0))*u.degree, frame='icrs')
    down1_left1 = SkyCoord(ra=(center_ra+((shift)*1.0/900.0/(np.cos(np.degrees(down1.dec)))))*u.hour, dec=down1.dec, frame='icrs')
    down1_right1 = SkyCoord(ra=(center_ra-((shift)*1.0/900.0/(np.cos(np.degrees(down1.dec)))))*u.hour, dec=down1.dec, frame='icrs')

    return(grid)

def spiral(X, Y):
    x = y = 0
    dx = 0
    dy = -1
    delta = np.zeros((max(X, Y)**2,2))
    for i in range(max(X, Y)**2):
        if (-X/2 < x <= X/2) and (-Y/2 < y <= Y/2):
            delta[i] = np.array([x,y])
        if x == y or (x < 0 and x == -y) or (x > 0 and x == 1-y):
            dx, dy = -dy, dx
        x, y = x+dx, y+dy
    return delta