from scipy import *
from scipy import optimize, linalg, integrate, interpolate
import pylab as pl

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

#def vermodo

def form(building):
	
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
	m = zeros((nfloors,nfloors))

	for i in range(nfloors):
		if i == (nfloors - 1):
			k[i,i] = kcol[i]
		else:
			k[i,i] = kcol[i] + kcol[i+1]
			k[i,i+1] = -kcol[i+1]
			k[i+1,i] = -kcol[i+1]
			
		m[i,i] = m0[i]
	
	
	# Find modes
	D,PHI = linalg.eig(k,m);
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
	building['m'] = m
	building['nfloors'] = nfloors
	building['rsis'] = ones((nfloors,1)) #Seismic input matrix
	
	return building
		
def plotmode(building,mode):
	th = arange(0., 2*pi, pi/20.)
	z = building['h'].cumsum()
	esp = building['esp']
	idx = building['modeorder']
	
	pl.figure()
	pl.plot(0*z,z,'--k')
	for i in range(building['nfloors']):
		rx = esp[i]*cos(th)
		ry = esp[i]*sin(th)
		factor = 0.1*z.max()/abs(building['phi'][idx[mode-1],:]).max()
		pl.plot(rx + factor*building['phi'][idx[mode-1],i], ry + z[i],'b')
	pl.xlabel('Disp [m]')
	pl.ylabel('z [m]')
	pl.title('Mode Num. {0:.0f}, T = {1:.2f} [s]'.format(mode,building['T'][idx[mode-1]]))
	pl.axis('equal')
	pl.show()
	
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

def freqresp(building,type,f0,f1,nfreq,plottype = 'plot'):
	
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
	
	pl.figure()
	ax = pl.subplot(1,1,1)
	for i in arange(building['nfloors']):
		
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
	pl.show()
