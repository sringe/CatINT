from transport.transport import Transport
from transport.calculator import Calculator
from transport.plot import Plot

tp=Transport()
tp.set_calculator('FTCS')
ca=Calculator(transport=tp)
pt=Plot(transport=tp)
pt.plot()
