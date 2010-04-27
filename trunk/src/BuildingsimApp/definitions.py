from scipy import *
from scipy import optimize, linalg, integrate, interpolate
import pylab as pl
import time

def findperiod(k, beta, Dx, Dy, ex, ey, m0):
    kx = k
    ky = kx/beta

    J0 = (Dx**2 + Dy**2)/12
    M0 = matrix([[m0, 0, -m0*ey],[0, m0, m0*ex],[-m0*ey, m0*ex, m0*J0]])
    kth = Dx*ky + Dy*kx
    K0 = diag(array([kx, ky, kth]))

    D, PHI = linalg.eig(K0,M0)
    T = 2*pi/sqrt(real(D))
    return T.max()

def objfun(k, beta, Dx, Dy, ex, ey, m0,Tobj):
    return Tobj - findperiod(k, beta, Dx, Dy, ex, ey, m0)

def msize(*args):#MATLAB-like size function
    Nargs = len(args)
    sz = array(args[0].shape)
    #return Nargs
    if Nargs == 1:
        return sz
        
    if Nargs > 1:
        return sz[args[1]]

def blkdiag(*args):
    Nargs = len(args)
    A = args[0]
    for i in range(Nargs-1):
        b = zeros((msize(A,0),msize(args[i+1],1)))
        A = vstack([hstack([A,b]),hstack([b.transpose(),args[i+1]])])
    return A
def form(building,geom_eff=0):
	
	g = 9.806
	
	E = building['E']
	h = building['h']
	I = building['I']
	dx = building['dx']
	dy = building['dy']
	gamma = building['gamma']
	esp = building['esp']
	xsi = building['xsi']
	
	#Calculate column stiffness
	kcol = 4*12*E*I/(h**3)
	nfloors = h.shape[0]
	m0 = dx*dy*esp*gamma/g     	# Main system mass tonf*s^2/m

	k = zeros((nfloors,nfloors))
	kg = zeros((nfloors,nfloors))
	m = zeros((nfloors,nfloors))
	
	if geom_eff == 1:
		P = -m0[::-1].cumsum()[::-1]*g

	for i in range(nfloors):
		if i == (nfloors - 1):
			k[i,i] = kcol[i]
		
		else:
			k[i,i] = kcol[i] + kcol[i+1]
			k[i,i+1] = -kcol[i+1]
			k[i+1,i] = -kcol[i+1]
		
		if geom_eff == 1:
			if i == (nfloors - 1):
				kg[i,i] = P[i]/h[i]
			else:
				kg[i,i] = P[i]/h[i] + P[i+1]/h[i+1]
				kg[i,i+1] = -P[i+1]/h[i+1]
				kg[i+1,i] = -P[i+1]/h[i+1]
				
		m[i,i] = m0[i]
	
	
	# Find modes
	D,PHI = linalg.eig(k+kg,m);
	w = sqrt(real(D));
	T = 2.*pi/w
	idx = argsort(T) #ascending period sort
	idx = idx[::-1]  #reverse order (to descending period)

	#Form xsi% modal classical damping matrix
	Mmodal = mat(dot(PHI.T,dot(m,PHI)))
	PHI = mat(PHI)
	M0 = mat(m)
	C0 = M0*PHI*(linalg.inv(Mmodal))*(2*xsi/100.*diag(w)*Mmodal)*linalg.inv(Mmodal)*PHI.T*M0

	#Actualizar estructura
	building['T'] = T
	building['c'] = C0
	building['phi'] = PHI
	building['modeorder'] = idx
	building['k'] = k	
	building['kg'] = kg
	building['m'] = m
	building['nfloors'] = nfloors
	building['rsis'] = ones((nfloors,1)) #Seismic input matrix
	if geom_eff==1:
		building['P'] = P
	
	return building
		
def plotmode(building,mode,anim=0,Nframe = 100,fps = 30.):
	z = building['h'].cumsum()
	esp = building['esp']
	idx = building['modeorder']
	
	factor = 0.1*z.max()/abs(building['phi'][idx[mode-1],:]).max()
	
	pl.ion()
	pl.figure()
	iter = 0
	fac = factor*cos(2*pi*iter/(2*fps))
	han = plotdef(building,squeeze(array(fac*building['phi'][idx[mode-1],:])))
	pl.xlabel('Disp [m]')
	pl.ylabel('z [m]')
	pl.title('Mode Num. {0:.0f}, T = {1:.2f} [s]'.format(mode,building['T'][idx[mode-1]]))
	pl.axis('equal')
	
	while iter < Nframe:
		if anim == 0:
			break
		time.sleep(1./fps)
		fac = factor*cos(2*pi*iter/(2*fps))
		plotdef(building, squeeze(array(fac*building['phi'][idx[mode-1],:])), update=1, handles=han)
		iter += 1
			

def animdef(building,u,Nframe = 1,factor=1,dt=1./30,fps = 30.):
	z = building['h'].cumsum()
	esp = building['esp']
	idx = building['modeorder']
	
	pl.ion()
	pl.figure()
	han = plotdef(building,squeeze(array(factor*u[0,:])))
	pl.xlabel('Disp [m]')
	pl.ylabel('z [m]')
	#pl.title('t = {0:.2f} [s]'.format(t[i]))
	pl.axis('equal')
	iter = 1
	fskip = ceil(dt*fps)
	if fskip < 1:
		fskip = 1
	if dt*fps < 1.:
		fskip = 1.
		fps = 1./dt
		
	
	while iter < min(Nframe*fskip,Nframe):
		time.sleep(1./fps)
		#ti = time.time()
		pl.title('t = {0:05.2f} s'.format(dt*iter))
		plotdef(building, squeeze(array(factor*u[iter,:])), update=1, handles=han)
		iter += fskip
		#tf = time.time()
		#if (tf-ti) > 1/fps:
			#fskip = floor((tf-ti)*fps)


def plotdef(building,u,update=0,handles={}):
	Npts = 30
	esp = building['esp']
	dx = building['dx']
	z = building['h'].cumsum()
	
	#Interpolation functions
	s = arange(0.,1.,1./Npts)
	N1 = 1.-3.*s**2.+2*s**3.
	N2 = 3.*s**2.-2.*s**3.
	flr = {}
	col = {}
	
	if update != 0:
		flr = handles['flrs']
		col = handles['cols']
	xx = array([-dx/2, dx/2 , dx/2, -dx/2])
	for i in range(building['nfloors']):		
		yy = array([-esp[i]/2, -esp[i]/2 , esp[i]/2, esp[i]/2])
		if update == 0:
			flr[i], = pl.fill(xx + u[i], z[i] + yy, facecolor = '#ADADAD', alpha=0.7, edgecolor='k')
		else:
			ux = xx + u[i]
			uy = z[i] + yy
			xy = flr[i].get_xy()
			xy[:,0] = ux[[0,1,2,3,0]]
			flr[i].set_xy(xy)
			
		if i == 0:
			z1 = 0
			z2 = z[i]
			u1 = 0
			u2 = u[i]
		else:
			z1 = z[i-1]
			z2 = z[i]
			u1 = u[i-1]
			u2 = u[i]
		if update == 0:
			col[2*i-1], = pl.plot(-dx/2 + N1*u1 + N2*u2, z1 + (z2-z1)*(s),'b')
			col[2*i],   = pl.plot(dx/2 + N1*u1 + N2*u2, z1 + (z2-z1)*(s),'b')
		else:
			col[2*i-1].set_xdata(-dx/2 + N1*u1 + N2*u2)
			col[2*i].set_xdata(dx/2 + N1*u1 + N2*u2)
	pl.draw()
	handles = dict(flrs=flr,cols=col)
	return handles
		

def response(building, input):
	#State-space representation

	A = mat(hstack([-linalg.solve(building['m'],building['k']),-linalg.solve(building['m'],building['c'])]) )
	b = -building['rsis']
	Ndof = building['nfloors']
	z0 = zeros((2*Ndof),'d') #Starting condition
	ug_f = interpolate.interp1d(input['t'],input['ug'],fill_value=0, bounds_error=0)

	#argtuple = (A,b,ug_f)
	
	#function for state-space representation of system
	
	def sssys(t, y):
		N = len(y)
		dzdt = 0.*y
		dzdt[:(N/2)] = y[(N/2):]
		dzdt[(N/2):] = dot(y,A.T) + ug_f(t)*b[:,0].T

		return squeeze(dzdt).real

	#Integration options
	odesolver = integrate.ode(sssys)
	odesolver.set_integrator('vode', method='bdf', order=5)
	odesolver.set_initial_value(z0)

	tf = input['t'].max()
	tout = input['t']
	dt = tout[2] - tout[1]
	zout = zeros((len(input['t']),msize(A,1)))
	i = 0
	while odesolver.successful() and odesolver.t < tf:
		odesolver.integrate(odesolver.t + dt)
		tout[i] = odesolver.t
		zout[i,:] = odesolver.y
		i += 1
	
	sol = {}
	sol['dis'] = zout[:,:Ndof]
	sol['vel'] = zout[:,Ndof:]
	sol['acc'] = (A*matrix(zout.T) + b*matrix(ug_f(tout))).T
	sol['t'] = tout
	
	return sol

def floorresp(building,input,sol,floor):

	pl.figure()
	pl.subplot(4,1,1)
	pl.plot(input['t'],input['ug'])
	pl.ylabel('Ug')
	pl.title('Response at floor Num. {0: .0f}'.format(floor))

	pl.subplot(4,1,2)
	pl.plot(sol['t'],sol['dis'][:,floor-1])
	pl.ylabel('D [m]')
	pl.subplot(4,1,3)
	pl.plot(sol['t'],sol['vel'][:,floor-1])
	pl.ylabel('V [m/s]')
	pl.subplot(4,1,4)
	pl.plot(sol['t'],sol['acc'][:,floor-1])
	pl.ylabel('A [m/**s]')
	pl.xlabel('t [sec]')

	pl.show()
	
def envresp(building,input,sol,type='acc'):

	z = building['h'].cumsum()
	pl.figure()
	#factor = 0.1*z.max()/abs(building['acc']).max()
	pl.plot(0*z,z,'k--')
	for n in arange(building['nfloors']):
		xord = sol[type][:,n]
		yord = z[n] + 0*xord
		pl.plot(xord,yord,'r')
	
	
	xmax = sol[type].max(0)
	xmin = sol[type].min(0)
	pl.plot(squeeze(array(xmax)),squeeze(z),'b')
	pl.plot(squeeze(array(xmin)),squeeze(z),'b')
	
	pl.show()

def freqresp(building,type=0,f0=0.,f1=10.,nfreq=200.,plottype = 'plot', plotthese = -1, axhan=-1):
	
	M = matrix(building['m'])
	C = building['c']
	K = building['k']
	r = matrix(building['rsis'])
	
	omegas = arange(2*pi*f0,2*pi*f1,2*pi*(f1-f0)/nfreq)
	U = zeros((building['nfloors'],nfreq))
	
	for i in arange(nfreq):
		w = omegas[i]
		U[:,i] = abs((1j*w)**type * linalg.solve(-w**2*M + (1j*w)*C + K, -M*r)).T
	
	#freq = {}
	#freq['f'] = omegas/(2*pi)
	#freq['U'] = U
	#return freq
	
	if plotthese == -1:
		plotthese = arange(building['nfloors'])
	
	if axhan==-1:
		pl.figure()
		ax = pl.subplot(1,1,1)
	else:
		ax = axhan
	
	for i in plotthese:
		
		if plottype.lower() == 'loglog':
			ax.loglog(squeeze(omegas)/(2*pi), squeeze(U[i,:]), label='Piso {0:.0f}'.format(i))
		elif plottype.lower() == 'semilox':
			ax.semilogx(squeeze(omegas)/(2*pi), squeeze(U[i,:]), label='Piso {0:.0f}'.format(i))
		elif plottype.lower() == 'semilogy':
			ax.semilogy(squeeze(omegas)/(2*pi), squeeze(U[i,:]), label='Piso {0:.0f}'.format(i))
		else:
			ax.plot(squeeze(omegas)/(2*pi), squeeze(U[i,:]), label='Piso {0:.0f}'.format(i))
			
	handles, labels = ax.get_legend_handles_labels()
	ax.legend(handles, labels)
	pl.xlabel('f [Hz]')
	pl.title('Funcion de respuesta en frecuencia')
	if type == 0:
		pl.ylabel('Desplazamientos')
	elif type == 1:
		pl.ylabel('Velocidades')
	elif type == 2:
		pl.ylabel('Aceleraciones')
	
def compare_buildings(build1,build2,type=0,f0=0.,f1=20.,nfreq=200.,plottype = 'plot'):
	Nf1 = build1['nfloors']
	Nf2 = build2['nfloors']
	print ''
	print '----------------------------------------'
	print 'Building Comparison'
	print '----------------------------------------'
	print 'Building 1 Name: '+build1['name']
	print 'Building 2 Name: '+build2['name']
	print '----------------------------------------'
	print '           Building 1     Building 2    '
	print 'Nfloors ={0:15.0f}{1:15.0f}'.format(Nf1,Nf2)
	print 'Mass    ={0:15.5f}{1:15.5f}'.format(build1['m'].sum(),build2['m'].sum())
	print 'Weight  ={0:15.5f}{1:15.5f}'.format(build1['m'].sum()*build1['g'],build2['m'].sum()*build1['g'])
	print ''
	print 'Periods [s]'

	for i in range(max(Nf1,Nf2)):
		if i > Nf1:
			print 'Mode{0:3.0f} =               {1:15.7f}'.format(i+1,build2['T'][build2['modeorder'][i]])
		elif i > Nf2:
			print 'Mode{0:3.0f} ={1:15.7f}               '.format(i+1,build1['T'][build1['modeorder'][i]])
		else:
			print 'Mode{0:3.0f} ={1:15.7f}{2:15.7f}'.format(i+1,build1['T'][build1['modeorder'][i]],build2['T'][build2['modeorder'][i]])
			
	print ''
	print 'Frequencies [Hz]'

	for i in range(max(Nf1,Nf2)):
		if i > Nf1:
			print 'Mode{0:3.0f} =               {1:15.7f}'.format(i+1,1/build2['T'][build2['modeorder'][i]])
		elif i > Nf2:
			print 'Mode{0:3.0f} ={1:15.7f}               '.format(i+1,1/build1['T'][build1['modeorder'][i]])
		else:
			print 'Mode{0:3.0f} ={1:15.7f}{2:15.7f}'.format(i+1,1/build1['T'][build1['modeorder'][i]],1/build2['T'][build2['modeorder'][i]])
	
	fig = pl.figure()	
	ax1 = pl.subplot(1,2,1)
	freqresp(build1, axhan = ax1, type = type, f0 = f0, f1 = f1, nfreq = nfreq, plottype = plottype)
	pl.title(build1['name'])
	
	ax2 = pl.subplot(1,2,2)
	freqresp(build2, axhan = ax2, type = type, f0 = f0, f1 = f1, nfreq = nfreq, plottype = plottype)
	pl.ylabel('')
	pl.title(build2['name'])
	
	fmax = max(ax1.axis()[1],ax2.axis()[1])
	fmin = min(ax1.axis()[0],ax2.axis()[0])
	ymax = max(ax1.axis()[3],ax2.axis()[3])
	ymin = min(ax1.axis()[2],ax2.axis()[2])
	
	ax1.axis(array([fmin,fmax,ymin,ymax]))
	ax2.axis(array([fmin,fmax,ymin,ymax]))
	
	
	handles = {}
	handles['fig'] = fig
	handles['ax1'] = ax1
	handles['ax2'] = ax2
	
	return handles

