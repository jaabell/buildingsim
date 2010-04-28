#!/usr/bin/python

from scipy import *
from scipy import optimize, linalg, integrate, interpolate
import time
from definitions import *


#columnas rectangulares
c1 = 30./100.                       # [m]
c2 = 30./100.                       # [m]
I0 = c1*(c2**3)/12.                 # [m**4]

# Global system parameters
building = {}
building['dx'] = 10.                                # [m]
building['dy'] =  6.                                # [m]
building['g'] = 9.806                               # Gravitational Constant [m/s^2]
building['xsi'] = 5.                                # Modal damping  [%]
building['h'] = array([2.7, 2.7, 2.7, 2.7])         # Altura entrepiso [m]
building['E'] = array([230., 230., 230., 230.])     # Modulo elasticidad piso [tonf/cm**2]
building['I'] = array([I0, I0, I0, I0]) * (100**2)  # Momento inercia columnas [m**4]
building['gamma'] = array([2.5, 2.5, 2.5, 2.5])     # Peso unitario losas [tonf/m**3]
building['esp'] = array([0.15, 0.15, 0.10, 0.10])   # espesor losas [m]
building['cadd'] = array([])

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



#han = compare_buildings(building1,building2,plottype='semilogy',type=2)
#pl.ion()
#pl.show()
sol = response(building,input)


#animdef(building,sol['dis'],Nframe=nt,factor=50,dt=dt,fps=2/dt)
#plotmode(building,1,anim=1)

#for i in arange(4):
    #plotmode(building,i+1)
#pl.show()
#plotdef(building,array([1, 2, 3, 4]))
#pl.show()


###Plot outputs
#
han = envresp(building,input,sol,type='dri')
pl.ion()
pl.show()
#envresp(building,input,sol,'dis')

#freqresp(building,0.,0.,20.,100.,'semilogx')
#freqresp(building,2.,0.,20.,100.,'semilogx')

