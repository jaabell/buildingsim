from scipy import *
#Harmonic or Chirp-Like
def harmonic(t,a,f=1.,b=0,K=0.,tmin=0.,tmax=Inf,):
	y = (a*sin(2.*pi*(f*t + K*t**2/2)) + b*cos(f*t + K*t**2/2))*(t>=tmin)*(t<=tmax)
	return y

#Rectangular
def rectangular(t,a,t0=0.,f=1.,tmin=0.,tmax=Inf):
	y = a*sign(sin(2*pi*f*(t-t0)))*(t>=tmin)*(t<=tmax)
	return y
	
#Triangular
#def triangular(t,a,t0=0.,f=1.,tmin=0.,tmax=Inf):
	#t_ = 
	#y = a*sign(sin(2*pi*f*(t-t0)))*(t>=tmin)*(t<=tmax)
	#return y	
	
#def stepsine(t,a,durs,freqs)
	#N = 
	#return y
#def 
