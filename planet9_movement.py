# -*- coding: utf-8 -*-
"""
Created on Tue Nov  1 15:07:13 2016

@author: Liam

02Nov2016 (jswift): Eliminated "math" module. 
                    Created default keywords instead of defining values internally
                    Adjusted number of days in a year to 365.2422
                    
14Nov2016 (Liam):   Created "find cadence" function
                    fixed error that said answer was in arcsecconds/min (its in arcessonds / hour)                
"""
#Planet 9 Location Finder

import numpy as np

#given summer solsace 2016 is 194 days from January 1 2017

def planet_speed(t,ra=90,dec=20,a=800):
    """
    Return proper motion of planet for given input parameters
    t: day of year
    ra: Right Ascencion in degrees
    dec: Declination in degrees
    a: Semi-major axis of planet (not accounting for Earth distance).
    """
    # Wrangle input variables
    RA = np.radians(np.float(ra))
    Dec = np.radians(np.float(dec))
    Dis = np.float(a)

    #calculates width and height of planet 9's movement.
    #Width and height are the distance from the ceter to poles, not absolute width / height
    width = (3600.0*180.0)/(Dis*np.pi)
    height = np.sin(Dec)*width

    #takes the derivative of x and y movement (parabolic equations)
    x_deriv = width*2*np.pi/365.2422*np.cos((t+195)*2*np.pi/365.2422)
    y_deriv = -height*2*np.pi/365.2422*np.sin((t+195)*2*np.pi/365.2422)

    #Calculates speed of planet 9 in Arcsecconds / hour
    speed = (((x_deriv)**2+(y_deriv)**(2))**(.500))/24
    return speed

def find_cadence(t,ra=90,dec=20,a=800,m=10):
    
    # Wrangle input variables
    RA = np.radians(np.float(ra))
    Dec = np.radians(np.float(dec))
    Dis = np.float(a)
    movement = np.float(m)

    #calculates width and height of planet 9's movement.
    #Width and height are the distance from the ceter to poles, not absolute width / height
    width = (3600.0*180.0)/(Dis*np.pi)
    height = np.sin(Dec)*width

    #takes the derivative of x and y movement (parabolic equations)
    x_deriv = width*2*np.pi/365.2422*np.cos((t+195)*2*np.pi/365.2422)
    y_deriv = -height*2*np.pi/365.2422*np.sin((t+195)*2*np.pi/365.2422)

    #Calculates speed of planet 9 in Arcsecconds / hour
    speed = (((x_deriv)**2+(y_deriv)**(2))**(.500))/24
    
    
    #Defines cadence in hours
    cadence = 10 / (speed*24)
    
    return cadence
    
    
    
    
    
    
    
    