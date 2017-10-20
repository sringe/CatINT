from transport.transport import Transport
from transport.calculator import Calculator 
from transport.plot import Plot
#t=Transport(integrator='Crank-Nicolson')


#odeint
#odeint--LF
#Crank-Nicolson
#Crank-Nicolson--LF
#FTCS
#FTCS--LF



pb_bound={
        'potential': {'wall':'zeta','bulk':0.0}}
#        'gradient': {'bulk':0.0}}
#pb_bound={
#        'potential': {'wall':'zeta','bulk':0.0}}

tp=Transport(pb_bound=pb_bound,nx=1000,system={'vzeta':-0.013}) #,nx=500)
tp.set_calculator('odeint')
tp.set_calculator('Crank-Nicolson--LF')
tp.set_initial_concentrations('Gouy-Chapman',vzeta=-0.012)

c=Calculator(transport=tp,tau_jacobi=1e-7,ntout=5,dt=1e-8,tmax=1e-5) #,scale_pb_grid='linear')
#scale_pb_grid
cout=c.run() #dt=1e-10,tmax=1e-5)

p=Plot(transport=tp)
p.plot(cout)
