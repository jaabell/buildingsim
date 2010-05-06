#!/usr/bin/python

from scipy import *
from scipy import optimize, linalg, integrate, interpolate
import time
from definitions import *


#columnas rectangulares
building = load_building('ejemplo8pisos.bsim')

g = 9.806

#Definicion input
a0 = 0.4*g
dt = 0.01                           # timestep
tmax = 30
tend = 20

fmin = 0.
fmax = 5.
timespan = 25.
f = 0.8
K = (fmax - fmin)/timespan

input = {}
t = arange(0., tmax , dt)    #Define time interval
nt = t.shape[0]
input['t'] = t
input['ug'] = a0*sign(sin(2.*pi*f*t))*(t<tend)                  #Rect
#input['ug'] = a0*sin(2.*pi*f*t)*(t<tend)                       #Sine
#input['ug'] = a0*sin(2.*pi*(f*t + K*t**2/2))*(t<tend)          #chirp

#d0 = 0.05
#sc = sin(2.*pi*(f*t + K*t**2/2))
#cc = cos(2.*pi*(f*t + K*t**2/2))
#input['ug'] = d0*(-sc*(2*pi*f + 2*pi*K*t)**2 + cc*2*pi*K)*(t<tend)             #test

building1 = building.copy()
building2 = building.copy()

building['name'] = 'Sin Efectos 2do. Orden'
building1['name'] = 'Sin Efectos 2do. Orden'
building2['name'] = 'Con Efectos 2do. Orden'
building2['gamma'] = building2['gamma']*10

#Initialize building and calculate response
building = form(building)
#building1 = form(building1,geom_eff=0)
#building2 = form(building2,geom_eff=1)

f0 = 0
f1 = 50
nfreq = 200
type = 0.

M = matrix(building['m'])
C = building['c']
K = building['k']
r = matrix(building['rsis'])
omegas = arange(2*pi*f0,2*pi*f1,2*pi*(f1-f0)/nfreq)
U = pl.complex128(zeros((building['nfloors'],nfreq)))

for i in arange(nfreq):
    w = omegas[i]
    U[:,i] = ((1j*w)**type * linalg.solve(-w**2*M + (1j*w)*C + K,-M*r)).T


#han = compare_buildings(building1,building2,plottype='semilogy',type=2)
#pl.ion()
#pl.show()
#sol = response(building,input)


#animdef(building,sol['dis'],Nframe=nt,factor=50,dt=dt,fps=2/dt)
#plotmode(building,1,anim=1)

#for i in arange(4):
    #plotmode(building,i+1)
#pl.show()
#plotdef(building,array([1, 2, 3, 4]))
#pl.show()


###Plot outputs
#
#han = envresp(building,input,sol,type='dri')
#pl.ion()
#pl.show()
#envresp(building,input,sol,'dis')

#freqresp(building,0.,0.,20.,100.,'semilogx')
#freqresp(building,2.,0.,20.,100.,'semilogx')

