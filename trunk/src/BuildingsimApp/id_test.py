#!/usr/bin/python

from scipy import *
from scipy import optimize, linalg, integrate, interpolate, fft
import time
import numpy as np
from definitions import *
from senales import *

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
dt = 0.01
tmax = 60

fmin = 0.
fmax = 20.
timespan = 30.
K = (fmax - fmin)/timespan

input = {}
t = arange(0., tmax , dt)    #Define time interval

ug = harmonic(t,a0,f=fmin,K = K, tmax=timespan)
input['t'] = t
input['ug'] = ug

building = form(building)
sol = response(building,input)

pl.ion()
#floorresp(building,input,sol,4)

u1 = sol['acc'][:,0]
u2 = sol['acc'][:,1]
u3 = sol['acc'][:,2]
u4 = sol['acc'][:,3]

w, U = transfun(building,type=2,f0=0,f1=20,nfreq=200.)
Ug = array(np.fft.fft(ug)).squeeze()

Nt = ug.size

U1 = array(fft(u1)).squeeze()
U2 = array(fft(u2)).squeeze()
U3 = array(fft(u3)).squeeze()
U4 = array(fft(u4)).squeeze()

Ue1 = abs(U1/Ug)
Ue2 = abs(U2/Ug)
Ue3 = abs(U3/Ug)
Ue4 = abs(U4/Ug)

Ureal1 = array(abs(U[0,:])).squeeze()
Ureal2 = array(abs(U[1,:])).squeeze()
Ureal3 = array(abs(U[2,:])).squeeze()
Ureal4 = array(abs(U[3,:])).squeeze()

wfft = np.fft.fftfreq(Nt,dt)

#pl.figure()
#pl.plot(np.fft.fftshift(wfft)/2./pi,abs(np.fft.fftshift(Ug)))


pl.figure()
pl.plot(np.fft.fftshift(wfft)/2/pi,np.fft.fftshift(Ue4))
pl.plot(w/2/pi,Ureal4)
pl.xlim([0,min(wfft.max(),w.max())/2/pi])


freqresp(building,type=2,f0=0,f1=20,nfreq=200.)
