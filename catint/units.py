import sys

#### some constants
unit_R=8.3144598
#kJ=kJ/mol
#kcal=kcal/mol
#J=Joule
#eV=electron volts
#AA=Angstrom
#Bohr=bohr radius
unit_e=1.6021766208e-19
eVToJ=unit_e #1.6021766208e-19
unit_eps0=8.854187817e-12
unit_NA=6.022140857e23
calToJ=4.182
unit_F=96485.33289
JTokcal = 1e-3/calToJ*unit_NA
BohrToAA = 0.52917721
HaToeV = 27.21138602
eVTokcal = 23.0609
eVTokcal2 = eVToJ*JTokcal
HaTokcal = 627.5095
HaToJ = 4.359744650e-18
unit_kB = 1.38064852e-23
unit_T = 298.14
unit_h = 6.626070040e-34
unit_c = 299792458
Rydberg=0.5*HaToeV
##################################

def convert_unit(val,start,end):
    #kcal/mol to x
    if start=='kcal' and end=='eV':
        factor=1./eVTokcal
    elif start=='kcal' and end=='kJ':
        factor=calToJ
    elif start=='kcal' and end=='kT':
        factor=1000./unit_R/unit_T*calToJ
    #eV to x
    elif start=='eV' and end=='kcal':
        factor=eVTokcal
    elif start=='eV' and end=='meV':
        factor=1000.
    elif start=='eV' and end=='Ha':
        factor=1./HaToeV
    elif start=='eV' and end=='kT':
        factor=eVTokcal*calToJ*1000./unit_R/unit_T
        #or: eVToJ/kB/T (but here 2 big variables)
    elif start=='eV' and end=='kJ':
        factor=eVTokcal*calToJ
        #or: eVToJ/1000.*NA
    #Ha to x
    elif start=='Ha' and end=='eV':
        factor=HaToeV
    elif start=='Ha' and end=='meV':
        factor=HaToeV*1000.
    elif start=='Ha' and end=='kJ':
        factor=HaTokcal*calToJ#HaToJ/1000.*NA
    elif start=='Ha' and end=='kcal':
        factor=HaTokcal
    elif start=='Ha' and end=='kT':
        factor=HaTokcal*calToJ*1000./unit_R/unit_T
        #or: HaToJ/kB/T
    #kJ/mol to x
    elif start=='kJ' and end=='kcal':
        factor=1./calToJ
    elif start=='kJ' and end=='kT':
        factor=1000./unit_R/unit_T
    #kT to x
    elif start=='kT':
        if end=='kJ':
            factor=1./1000.*unit_R*unit_T
        elif end=='kcal':
            factor=0.001*unit_R*unit_T/calToJ
        #or: JTokcal*kB*T
        elif end=='eV':
            factor=0.001*unit_R*unit_T/calToJ/eVTokcal
    #J to x
    elif start=='J':
        if end=='Ha':
            factor=1./HaToJ
        elif end=='eV':
            factor=1./eVToJ
        elif end=='kcal':
            factor=1./JTokcal
    elif start==end:
        factor=1.
    else:
        print('unexpected conversion units, check units module for available conversions')
        sys.exit()
    return val*factor
