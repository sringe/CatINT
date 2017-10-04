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

tp=Transport()
tp.set_calculator('odeint')

c=Calculator(transport=tp)
cout=c.run(dt=1e-9,tmax=1e-6)

p=Plot(transport=tp)
p.plot(cout)
