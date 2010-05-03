from scipy import *
#Harmonic or Chirp-Like
def harmonic(t,a,f=1.,b=0,K=0.,tmin=0.,tmax=Inf,):
    y = (a*sin(2.*pi*(f*t + K*t**2/2)) + b*cos(2.*pi*(f*t + K*t**2/2)))*(t>=tmin)*(t<=tmax)
    return y
#
def ensayo_real(t,intervalos,frecuencias,d,phi=0.):
    tpos = intervalos.cumsum()
    N = intervalos.size
    y = zeros(t.shape)
    for i in range(N):
        if i == 0:
            ti = 0
        else:
            ti = tpos[i-1]
        tf = tpos[i]
        f = frecuencias[i]
        y += -(2*pi*f)**2*d*sin(2*pi*f*t + phi)*( t >= ti)*(t<= tf)
        phi += 2*pi*f*intervalos[i]
        print 'ti = {0:5.1f}, tf = {1:5.1f}, f = {2:6.3f}'.format(ti,tf,f)
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
