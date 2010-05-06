# Contiene principales funciones de procesamiento y graficacion de 
# BuildingSimApp.
#
from scipy import *
from scipy import optimize, linalg, integrate, interpolate
import pylab as pl
import time
import numpy as np
import sys
import xml.etree.cElementTree as ET

#
#No se usa
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
#
#
#No se usa
def objfun(k, beta, Dx, Dy, ex, ey, m0,Tobj):
    return Tobj - findperiod(k, beta, Dx, Dy, ex, ey, m0)
#
#
#No se usa
def msize(*args):#MATLAB-like size function
    Nargs = len(args)
    sz = array(args[0].shape)
    #return Nargs
    if Nargs == 1:
        return sz

    if Nargs > 1:
        return sz[args[1]]
#
#
#De uso interno
# M = blkdiag(A,B,C,D,...)
# A,B,C,D,... son arreglos numpy. Esta funcion devuelve un arreglo M
# diagonal por bloques, donde los bloques son A,B,C,D,...
# La cantidad de arreglos no esta limitada. La matriz M creada es de 
# estructura llena (ceros explicitos).
def blkdiag(*args):
    
    Nargs = len(args)
    A = args[0]
    for i in range(Nargs-1):
        b = zeros((msize(A,0),msize(args[i+1],1)))
        A = vstack([hstack([A,b]),hstack([b.transpose(),args[i+1]])])
    return A
#
#
# FORM: Inicializa la estructura.
#
# Sintaxis: building = form(building)
# building: Estructura de datos (diccionario) que define al edificio.
#
# Building debe tener definidos los campos:
#   'E'     : Arreglo numpy con los modulos de elasticidad de cada piso
#   'I'     : Arreglo numpy con los momentos de inercia de las 
#             columnas de cada piso.
#   'dx'    : Dimension X de la planta del edificio.
#   'dy'    : Dimension Y de la planta del edificio.
#   'gamma' : Arreglo numpy con el peso unitario del material de la 
#             losa de cada piso.
#   'esp'   : Arreglo numpy con el espesor de la losa en cada piso.
#   'xsi'   : Amortiguamiento. (0 <= xsi < 1)
#   'name'  : String con el nombre del edificio
#
# Todos los valores son positivos.
#
# Form devuelve el mismo building con mas campos definidos (matrices 
# del sistema, etc.)
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
    m0 = dx*dy*esp*gamma/g      # Main system mass tonf*s^2/m

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
#
#
# PLOTMODE: Graficar forma modal
#
# Sintaxis:
# handle = plotmode(building,mode,anim=0,Nframe = 100,fps = 30.)
#
#   building: Estructura de datos (diccionario) que define al edificio.
#   mode    : (int) > 0 que indica cual modo graficar (por defecto = 1)
#   anim    : (bool) indica si animar el modo o no (por defecto no)
#   Nframe,fps : Parametros para la animacion. Numero de frames 
#                totales que animar a fps objetivo.
#
# Se puede ejecutar solo despues de haber aplicado  
#  building =form(building)
#
# Devuelve la handle a la figura generada
def plotmode(building,mode=1,anim=0,Nframe = 100,fps = 30.):
    z = building['h'].cumsum()
    esp = building['esp']
    idx = building['modeorder']

    factor = 0.1*z.max()/abs(building['phi'][idx[mode-1],:]).max()

    handle = pl.figure()
    iter = 0
    fac = factor*cos(2*pi*iter/(2*fps))
    han = plotdef(building,squeeze(array(fac*building['phi'][:,idx[mode-1]])))
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

    return handle
#
#
#ANIMDEF: Animar la deformada de la estructura
#
#   building: Estructura de datos (diccionario) que define al edificio.
#   u       : Arreglo numpy (matriz) con la solucion del problema.
#
#   Si,
#   sol = response(building,input).
#   u = sol['dis']
def animdef(building,u,Nframe = 1,factor=1,dt=1./30,fps = 30.):

    #Calcular alturas y definir algunas variables para facilitar el 
    #codigo
    z = building['h'].cumsum()
    esp = building['esp']
    idx = building['modeorder']

    #Entrar modo interactivo e inicializar deformada en cero
    pl.ion()
    pl.figure()
    han = plotdef(building,squeeze(array(factor*u[0,:])))
    pl.xlabel('Disp [m]')
    pl.ylabel('z [m]')
    pl.axis('equal')
    
    #Calcular fps y fskip si hay inconsistencias
    iter = 1
    fskip = ceil(dt*fps)
    if fskip < 1:
        fskip = 1
    if dt*fps < 1.:
        fskip = 1.
        fps = 1./dt

    #Ciclo principal de animacion
    while iter < min(Nframe*fskip,Nframe):
        time.sleep(1./fps)
        pl.title('t = {0:05.2f} s'.format(dt*iter))
        plotdef(building, squeeze(array(factor*u[iter,:])), update=1, handles=han)
        iter += fskip
    
    #Salir de modo interactivo y cerrar ventana de animacion
    pl.ioff()
    #pl.close(all)
#
#
#PLOTDEF: Graficar la deformada estatica de la estructura (sin animar)
#
#   Sintaxis:
#   handles = plotdef(building,u,update=0,handles={},ax = -1):
#
#   handles:  (output si update==0) Es un diccionario con los handles a 
#             la figura y elementos creados.
#
#   building: Estructura de datos (diccionario) que define al edificio.
#   u       : Arreglo numpy (vector) con la deformada a graficar.
#   update  : (bool) que indica si alcualizar una deformada existente 
#              o dibujar una nueva. (por defecto 0)
#   handles : handles de una figura ya creada (sirve si update == 1)
#
#  Ejemplo
#
#  sol1 = response(building,input1)
#  sol2 = response(building,input1)
#
#  Seleccionar respuestas para el instante t[10]
#  u1 = sol1['dis'][10,:]
#  u2 = sol1['dis'][10,:]
#
#  handles = plotdef(building,u1)  (crea la figura usando u1 como 
#                                   solucion )
#  Actualizar la ventana con la solucion u2.
#  plotdef(building,u2,update = 1, handles = handles)
#            
def plotdef(building,u,update=0,handles={},ax = -1):
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
    #if ax != -1:
        #pl.axes(ax)
    #else:
        #ax = subplot(111)
        
    if update != 0:
        flr = handles['flrs']
        col = handles['cols']
    xx = array([-dx/2, dx/2 , dx/2, -dx/2])
    for i in range(building['nfloors']):
        yy = array([-esp[i]/2, -esp[i]/2 , esp[i]/2, esp[i]/2])
        if update == 0:
            flr[i], = pl.fill(xx + u[i], z[i] + yy, facecolor = 
            '#ADADAD', alpha=0.7, edgecolor='k')
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
            col[2*i],   = pl.plot(dx/2 + N1*u1 + N2*u2, z1 +
            (z2-z1)*(s),'b')
        else:
            col[2*i-1].set_xdata(-dx/2 + N1*u1 + N2*u2)
            col[2*i].set_xdata(dx/2 + N1*u1 + N2*u2)
    #if ax == -1:
        #pl.draw()#
    handles = dict(flrs=flr,cols=col)
    return handles
#
#
# RESPONSE: Calcular solucion en el tiempo al problema definido.
#
#   building: Estructura de datos (diccionario) que define al edificio.
#   input   : Estructura de datos (diccionario) que define el input.
#
#  input debe tener los siguientes campos:
#       't'  : Arreglo numpy (vector) con el tiempo
#       'ug' : Arreglo numpy (vector) con las aceleracciones 
#              correspondientes
#   ug[i] es la aceleraccion correspondiente al tiempo t[i]
#
#  La simulacion se ejecuta para todo los tiempo definidos
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
#
#
# FLOORRESP: Graficar solucion en el tiempo para un piso en particular
#
# Sintaxis:
#  h = floorresp(building,input,sol,floor)
#
#   building: Estructura de datos (diccionario) que define al edificio.
#   input   : Estructura de datos (diccionario) que define el input.
#   sol     : Estructura de datos (diccionario) que define la solucion 
#             (generada por sol = response(building,input).
#   floor   : (int)>=1 que define el piso donde se quiere graficar la 
#              respuesta.
# 
#   h       : (output) handle a la figura generada
def floorresp(building,input,sol,floor):

    h = pl.figure()
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

    return h
#
#
# ENVRESP: Graficar solucion envolvente para algun tipo de output
#
# Sintaxis:
#   handles = envresp(building,input,sol,type='acc')
#
#   building: Estructura de datos (diccionario) que define al edificio.
#   input   : Estructura de datos (diccionario) que define el input.
#   sol     : Estructura de datos (diccionario) que define la solucion 
#   type    : (string) Tipo de grafico a generar, puede ser uno de
#       type = 'acc': Graficar envolventes de aceleraciones de piso
#       type = 'vel': Graficar envolventes de velocidades de piso
#       type = 'dis': Graficar envolventes de desplazamientos de piso
#       type = 'dri': Graficar envolventes de drift de entrepiso (implementar)
#       type = 'cor': Graficar envolventes de corte de entrepiso (implementar)
#       type = 'mom': Graficar envolventes de momento de entrepiso (implementar)
def envresp(building,input,sol,type='acc'):
    #Evaluar alturas
    h = building['h']
    z = h.cumsum()
    
    #Inicializar figura
    pl.figure()
    
    #factor = 0.1*z.max()/abs(building['acc']).max()   #Factor de escala
    
    #Calcular respuestas maximas
    if (type == 'acc') | (type == 'dis') | (type == 'vel'):
        xmax = sol[type].max(0)
        xmin = sol[type].min(0)

    #Calcular PGAs del input
    amax = input['ug'].max()
    amin = input['ug'].min()

    #Ampliar vector de alturas
    znew = hstack((0.,z))
    z = znew

    if type == 'acc':
        xmax = hstack((amax,array(xmax).squeeze()))
        xmin = hstack((amin,array(xmin).squeeze()))
        xl = 'Aceleraciones'
    elif type == 'dis':
        xmax = hstack((0.,array(xmax).squeeze()))
        xmin = hstack((0.,array(xmin).squeeze()))
        xl = 'Desplazamientos'
    elif type == 'vel':
        xmax = hstack((0.,array(xmax).squeeze()))
        xmin = hstack((0.,array(xmin).squeeze()))
        xl = 'Velocidades'
    elif type == 'cor':
        drift = hstack((matrix(sol['dis'][:,0]).T, diff(sol['dis'],1) ))
        drmax = drift.max(0)
        drmin = drift.min(0)
        xmax = zeros(2*building['nfloors'])
        xmin = xmax.copy()
        znew = xmax.copy()
        for i in range(building['nfloors']):
            xmax[2*i] = drmax[0,i]/h[i]*100.
            xmax[2*i+1] = drmax[0,i]/h[i]*100.
            xmin[2*i] = drmin[0,i]/h[i]*100.
            xmin[2*i+1] = drmin[0,i]/h[i]*100.
            znew[2*i] = z[i]
            znew[2*i+1] = z[i+1]
            #xmax = hstack((0.,array(xmax).squeeze()))
            #xmin = hstack((0.,array(xmin).squeeze()))
        z = znew
        xl = '% Drift'
        xl = 'Corte'
    elif type == 'mom':
        #xmax = hstack((0.,array(xmax).squeeze()))
        #xmin = hstack((0.,array(xmin).squeeze()))
        xl = 'Momento'
    elif type == 'dri':
        drift = hstack((matrix(sol['dis'][:,0]).T, diff(sol['dis'],1) ))
        drmax = drift.max(0)
        drmin = drift.min(0)
        xmax = zeros(2*building['nfloors'])
        xmin = xmax.copy()
        znew = xmax.copy()
        for i in range(building['nfloors']):
            xmax[2*i] = drmax[0,i]/h[i]*100.
            xmax[2*i+1] = drmax[0,i]/h[i]*100.
            xmin[2*i] = drmin[0,i]/h[i]*100.
            xmin[2*i+1] = drmin[0,i]/h[i]*100.
            znew[2*i] = z[i]
            znew[2*i+1] = z[i+1]
            #xmax = hstack((0.,array(xmax).squeeze()))
            #xmin = hstack((0.,array(xmin).squeeze()))
        z = znew
        xl = '% Drift'
        
    
    #Panel izquierdo
    sp1 = pl.subplot(1,3,1)
    pl2 = pl.plot(squeeze(array(xmin)),squeeze(z),'b')
    pl.ylabel('Altura')
    pl.xlabel(xl)

    #Graficar edificio al centro
    sp2 = pl.subplot(1,3,2)
    plbuild = plotdef(building,0*z,ax=sp1)
    pl.title(building['name'])
        
    #Panel derecho
    sp3 = pl.subplot(1,3,3)
    pl3 = pl.plot(squeeze(array(xmax)),squeeze(z),'b')
    pl.xlabel(xl)

    #Ajustar ventanas
    x1 = sp1.get_position().xmax
    x2 = sp3.get_position().xmin
    w = x2 - x1
        
    sp2.set_position([x1, sp2.get_position().ymin, w, sp2.get_position().height])
    sp2.set_yticklabels('')
    sp2.set_xticklabels('')
    sp3.set_yticklabels('')
    sp1.set_ylim(sp2.get_ylim())
    sp3.set_ylim(sp2.get_ylim())
    xl = sp1.get_xlim()
    sp1.set_xlim([xl[0], 0])
    xl = sp3.get_xlim()
    sp3.set_xlim([0, xl[1]])
    
    #Devolver handles a las componentes de la figura
    handles = {}
    handles['plbuild'] = plbuild
    handles['pl2'] = pl2
    handles['pl3'] = pl3
    handles['sp1'] = sp1
    handles['sp2'] = sp2
    handles['sp3'] = sp3

    return handles
#
#
# FREQRESP: Graficar funcion de respuesta en frecuencia del sistema
#
# Sintaxis:
#   handle = freqresp(building,type=0,f0=0.,f1=10.,nfreq=200.,plottype = 'plot', plotthese = -1, axhan=-1)
#
#   building: Estructura de datos (diccionario) que define al edificio.
#   type    : (int) que define el tipo de output buscado. Puede ser:
#           type = 0: Output desplazamientos
#           type = 1: Output velocidades
#           type = 2: Output aceleraciones
#   f0      : Frecuencia inicial (minima) a graficar.
#   f1      : Frecuencia final (maxima) a graficar.
#   nfreq   : Numero de frecuencias a generar
#   plottype: String con el nombre de la funcion a usar para graficar. 
#       Por ejemplo:
#           plottype = 'plot'       :  Plot simple.
#           plottype = 'loglog'     :  Plot logaritmico doble.
#           plottype = 'semilogx'   :  Plot semilogaritmico en x.
#           plottype = 'semilogy'   :  Plot semilogaritmico en y.
#   plotthese: Arreglo numpy con (int)s que definen para que pisos se 
#              graficara.
#   axhan   : Handle a la figura en el que se desea graficar. (se usa 
#             al comparar funciones)
#
#  Devuelve handle a la figura generada.
def freqresp(building,type=0,f0=0.,f1=10.,nfreq=200.,plottype = 'plot', plotthese = -1, axhan=-1):

    omegas,Udum = transfun(building,type,f0,f1,nfreq)
    U = pl.absolute(Udum)

    if plotthese == -1:
        plotthese = arange(building['nfloors'])

    if axhan==-1:
        handle = pl.figure()
        ax = pl.subplot(1,1,1)
    else:
        ax = axhan
        handle = pl.gcf()

    for i in plotthese:
        if plottype.lower() == 'loglog':
            ax.loglog(squeeze(omegas)/(2*pi), squeeze(U[i,:]),
            label='Piso {0:.0f}'.format(i+1))
        elif plottype.lower() == 'semilogx':
            ax.semilogx(squeeze(omegas)/(2*pi), squeeze(U[i,:]),
            label='Piso {0:.0f}'.format(i+1))
        elif plottype.lower() == 'semilogy':
            ax.semilogy(squeeze(omegas)/(2*pi), squeeze(U[i,:]),
            label='Piso {0:.0f}'.format(i+1))
        else:
            ax.plot(squeeze(omegas)/(2*pi), squeeze(U[i,:]),
            label='Piso {0:.0f}'.format(i+1))

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
    
    return handle, U
#
#  TRANSFUN: Calcular funcion de transferencia.
#
# Misma entrada que freqresp, devuelve U.
def transfun(building,type=0,f0=0.,f1=10.,nfreq=200.):
    M = matrix(building['m'])
    C = building['c']
    K = building['k']
    r = matrix(building['rsis'])

    omegas = arange(2*pi*f0,2*pi*f1,2*pi*(f1-f0)/nfreq)
    U = pl.complex128(zeros((building['nfloors'],nfreq)))

    for i in arange(nfreq):
        w = omegas[i]
        U[:,i] = ((1j*w)**type * linalg.solve(-w**2*M + (1j*w)*C + K,-M*r)).T
    return omegas,U
#
#
# COMPARE_BUILDINGS: Comparar datos dinamicos y funciones de respuesta 
# en frecuencia para dos edificios.
#
# Sintaxis: 
#  compare_buildings(build1,build2,type=0,f0=0.,f1=20.,nfreq=200.,plottype = 'plot')
#
#   build1: Estructura de datos (diccionario) que define al edificio 1.
#   build2: Estructura de datos (diccionario) que define al edificio 2.
# 
# El resto de los parametros igual que FREQRESP. Importante que esta 
# funcion solo se puede usar despues de inicializar los edificios con:
#
#  build1 = form(build1)
#  build2 = form(build2)
#
# Ademas esta funcion escribe en la linea de comando algunos datos 
# para comparar ambos edificios. (implementar que la salida sea a un 
# .txt).
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
    sys.stdout.flush()

    for i in range(max(Nf1,Nf2)):
        if i+1 > Nf1:
            print 'Mode{0:3.0f} =               {1:15.7f}'.format(i+1,build2['T'][build2['modeorder'][i]])
        elif i+1 > Nf2:
            print 'Mode{0:3.0f} ={1:15.7f}               '.format(i+1,build1['T'][build1['modeorder'][i]])
        else:
            print 'Mode{0:3.0f} ={1:15.7f}{2:15.7f}'.format(i+1,build1['T'][build1['modeorder'][i]],build2['T'][build2['modeorder'][i]])

    print ''
    print 'Frequencies [Hz]'
    sys.stdout.flush()

    for i in range(max(Nf1,Nf2)):
        if i+1 > Nf1:
            print 'Mode{0:3.0f} =               {1:15.7f}'.format(i+1,1/build2['T'][build2['modeorder'][i]])
        elif i+1 > Nf2:
            print 'Mode{0:3.0f} ={1:15.7f}               '.format(i+1,1/build1['T'][build1['modeorder'][i]])
        else:
            print 'Mode{0:3.0f} ={1:15.7f}{2:15.7f}'.format(i+1,1/build1['T'][build1['modeorder'][i]],1/build2['T'][build2['modeorder'][i]])
            
    sys.stdout.flush()

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
#
# PLOTINPUT: Graficar un input
#  Recibe un input.
def plotinput(input):
    fig = pl.figure()
    fig.set_figwidth(11)
    fig.set_figheight(4)
    pl.plot(input['t'],input['ug'])
    pl.xlabel('t')
    pl.ylabel('ug')
    return fig

#
# Funciones de utilidad (nula)
def shiftplot(w,f):
    pl.plot(np.fft.fftshift(w),np.fft.fftshift(f))


def load_building(filename):
        tree = ET.parse(filename)
        root = tree.getroot()
        
        building = {}
        building['name'] = root.get("name")
        
        building['dx'] = float(root.get("dx"))
        building['dy'] =  float(root.get("dy"))
        building['g'] = 9.806				# Gravitational Constant [m/s^2]
        building['xsi'] = 5. 				# Modal damping  [%]        
                
        numfloors = sum(1 for floor in tree.getiterator("floor"))
        
        building['h'] = zeros((numfloors,))       # Altura entrepiso [m]
        building['E'] = zeros((numfloors,))       # Modulo elasticidad piso [tonf/cm**2]
        building['I'] = zeros((numfloors,))       # Momento inercia columnas [m**4]
        building['gamma'] = zeros((numfloors,))   # Peso unitario losas [tonf/m**3]
        building['esp'] = zeros((numfloors,))     # espesor losas [m]
        
        for i, floorNode in enumerate(tree.getiterator("floor")):
            floor = {"h":floorNode.get("h"), "E":floorNode.get("E"), "I":floorNode.get("I"), "gamma":floorNode.get("gamma"), "esp":floorNode.get("esp")}
            for key,value in floor.iteritems():
                building[key][i] = float(value) 
        
        return building
