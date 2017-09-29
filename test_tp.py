from transport import Transport

#t=Transport(integrator='Crank-Nicolson')
#t=Transport(integrator='FTCS-odeint')
t=Transport(integrator='FTCS')
t.run()
t.plot()
