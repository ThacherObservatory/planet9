# -*- coding: utf-8 -*-
"""
Created on Tue Jan 31 09:46:15 2017

@author: Liam
"""

import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord

def grid():
    
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