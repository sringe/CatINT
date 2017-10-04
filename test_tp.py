from transport import Transport

#t=Transport(integrator='Crank-Nicolson')

t=Transport(
    integrator='odeint', #,
    integrator_correction=False, #True,
    dt=1e-9,
    tmax=1e-6)

t.run()
t.plot()
