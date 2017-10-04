from transport import Transport

#t=Transport(integrator='Crank-Nicolson')
t=Transport(
    integrator='odeint',
    integrator_correction=False, #True,
    dt=1e-12,
    tmax=1e-9)
#t=Transport(integrator='FTCS')
t.run()
t.plot()
