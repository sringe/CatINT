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

#pb_bound={
#        'potential': {'wall':'zeta'},
#        'gradient': {'bulk':0.0}}
pb_bound={
        'potential': {'wall':'zeta','bulk':0.0}}

tp=Transport(pb_bound=pb_bound)
#tp.set_calculator('odeint')
tp.set_calculator('Crank-Nicolson--LF')
#tp.set_initial_concentrations('Gouy-Chapman')

c=Calculator(transport=tp,tau_jacobi=1e-6)
#scale_pb_grid
cout=c.run(dt=1e-12,tmax=1e-10)

p=Plot(transport=tp)
p.plot(cout)
