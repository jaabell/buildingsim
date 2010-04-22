#!/usr/bin/python


from scipy import *
from scipy import optimize, linalg, integrate, interpolate
import time
from definitions import *

def building_sim(building):
    #columnas rectangulares
    c1 = 30./100.						# [m]
    c2 = 30./100.						# [m]
    I0 = c1*(c2**3)/12.					# [m**4]

    # Global system parameters
    g = 9.806

    #Definicion input
    a0 = 0.4*g
    dt = 0.01               			# timestep
    tmax = 30
    tend = 20

    fmin = 0.
    fmax = 5.
    timespan = 25.
    f = fmin
    K = (fmax - fmin)/timespan

    input = {}
    t = arange(0., tmax , dt)    #Define time interval
    input['t'] = t
    input['ug'] = a0*sign(sin(2.*pi*f*t))*(t<tend)      			#Rect
    input['ug'] = a0*sin(2.*pi*f*t)*(t<tend) 						#Sine
    input['ug'] = a0*sin(2.*pi*(f*t + K*t**2/2))*(t<tend) 			#chirp

    d0 = 0.05
    sc = sin(2.*pi*(f*t + K*t**2/2))
    cc = cos(2.*pi*(f*t + K*t**2/2))
    input['ug'] = d0*(-sc*(2*pi*f + 2*pi*K*t)**2 + cc*2*pi*K)*(t<tend) 			#test


    #Initialize building and calculate response
    building = form(building)
    sol = response(building,input)


    ##Plot outputs
    floorresp(building,input,sol,4)
    envresp(building,input,sol)
    envresp(building,input,sol,'dis')
    freqresp(building,0.,0.,20.,100.,'semilogx')
    freqresp(building,2.,0.,20.,100.,'semilogx')

