# -*- coding: utf-8 -*-
"""
Created on Tue Nov  1 15:07:13 2016

@author: Liam
"""
#Planet 9 Location Finder

import numpy as np
import math


#given summer solsace 2016 is 194 days from January 1 2017

def planet_speed(t):
    #define Right Accension, Declination, Distance to planet 9 (au) and Days from January 1
    RA = np.radians(90.0)
    Dec = np.radians(20.0)
    Dis = 800

    #calculates width and height of planet 9's movement.
    #Width and height are the distance from the ceter to poles, not absolute width / height
    width = 1.00/Dis/math.pi*180*3600
    height = math.sin(Dec)*width

    #takes the derivative of x and y movement (parabolic equations)
    x_deriv = width*2*math.pi/365.5*math.cos((t+195)*2*math.pi/365.5)
    y_deriv = -height*2*math.pi/365.5*math.sin((t+195)*2*math.pi/365.5)

    #Calculates speed of planet 9 in Arcsecconds / min
    speed = (((x_deriv)**2+(y_deriv)**(2))**(.500))/24
    return speed