
# Definition of Kinetics 

## Fixed species fluxes 

We start by introducing into the definition of fluxes in CatINT as constants within descriptor space. We look at HER and CO2R to CO in alkaline conditions. At the electrode, we define the two reactions. The dictionary key must be the product name as also appearing in the species dictionary:

    electrode_reactions=
        {
        'H2': {'reaction': 'H2O + 2 e- -> H2 + 2 OH-'}
        'CO': {'reaction': 'H2O + CO2 + 2 e- -> CO + 2 OH-'} 
        }
        
All appearing species must be included in the species list (in spite of H$_2$O and e$^-$):

    species=
        {
    	'H2': { 'symbol': 'H_2',
		        'name': 'hydrogen',
		        'diffusion': 4.50e-009, #(m^2/s)
		        'bulk concentration': 1e-4}, #(mol/m^3)
        'OH-': {'symbol': 'OH^-',
		        'name': 'hydroxide',
	            'diffusion': '5.273e-9', #(m^2/s)
	            'bulk concentration': 1e-7}, #(mol/m^3)
	    'CO2': {'symbol': 'CO_2',
	            'name': 'carbon dioxide',
	            'diffusion': 1.91e-009, #(m^2/s)
	            'bulk concentration': 1e-3} #(mol/m^3)
	    }
We can now add fluxes. These can be given as production/consumption rate (flux) of the individual species:

    species['H2']['flux'] = 1e-4 #(mol/s/m^2)
    species['CO']['flux'] = 1e-5 #(mol/s/m^2)

Alternatively, we can provide them as current densities:

    species['H2']['flux'] = 1e-4 #(mol/s/m^2)
    species['CO']['flux'] = 1e-5 #(mol/s/m^2)
    
Not, however, that only one flux should be defined per equation. The fluxes of the remaining species will be calculated automatically from the reaction equation, for example:

$j_\mathrm{OH^-} = 2 \times j_\mathrm{H_2} + 2 \times j_\mathrm{CO}$
$j_\mathrm{CO_2} = -j_\mathrm{CO}$

Alternatively, fluxes can be also given as a current density, e.g. if experimental data is at hand that should be tested. This can be defined as:

    species['H2']['current density'] = 1e-4 #(C/s/m^2=A/m^2)
    species['CO']['current density'] = 1e-5 #(C/s/m^2=A/m^2)

Internally, these will be then converged into fluxes via:

$j_\mathrm{OH^-} = -2\times i_\mathrm{OH^-}$
$j_\mathrm{CO} = -2\times i_\mathrm{CO}$

Finally, fluxes can be also given as faradaic efficiencies (FE):

    species['H2']['FE'] = 0.6 #(%)
    species['CO']['FE'] = 0.2 #(%)

This, however, requires to also define the total current density in the system dictionary as:

    system['current density'] = 1e-4

If a total current density is given, a check will be performed if the partial current density add up to the total current density. If they do not, the program will stop with a warning. Please define then an unknown species with the missing current density.

## Rate-Equations

In COMSOL (and only here, not for any other of the FD solvers), we can define fluxes also as equations depending e.g. on the local concentrations of the species and the local potential. We can e.g. define the flux of H$_2$ in terms of a Butler-Volmer relation:

$j_\mathrm{H_2}=\rho\cdot\theta_\mathrm{H}\cdot A\cdot \exp(-[G_a^\mathrm{eq}+\alpha\cdot F\cdot \eta]/RT)$,

where $\rho$ is the active site density, $\theta_\mathrm{H}$ is the coverage of H, $G_a$ is the activation barrier of the rate-determining step and $\alpha$ the transfer coefficient. $\eta$ is the overpotential. Against a reference electrode which does not vary with pH, this is just given as the difference between potential at the reaction plane and potential at the electrode:

$\eta = (\Phi^\mathrm{M}-\Phi^\ddagger)-(\Phi^\mathrm{M,eq}-\underbrace{\Phi^{\ddagger,\mathrm{eq}}}_{\approx 0})=(\Phi^\mathrm{M}-\Phi^\ddagger)-\Phi^\mathrm{M,eq}$.

We now define some COMSOL variables and parameters:

    comsol_params['Ga']=[str(-0.3*unit_F)+'[J/mol]','Adsorption barrier H']
    comsol_params['alpha']=['0.5','Transfer Coefficient']
    comsol_params['A']=['1.e13[1/s]','Exponential prefactor']
    comsol_params['rho']=['1e-05[mol/m^2]','Density of Active Sites']
    comsol_params['theta_max']=['0.4','Maximum Coverage']
    comsol_params['Lmol']=['1[l/mol]','unit conversion factor']
    comsol_params['Kads']=['1e-4','Adsorption equilibrium constant']
    comsol_params['phiEq']=['-0.1[V]','equilibrium potential']

What is missing now is to define the coverage of H. We can assume a Langmuir isotherm with a maximum coverage of $\theta_\mathrm{H}^\mathrm{max}$:

$\theta_\mathrm{H}=\frac{\sqrt{K_\mathrm{ads}\cdot a_\mathrm{H_2}}}{1+\sqrt{a_\mathrm{H_2}K_\mathrm{ads}}}\theta_\mathrm{H}^\mathrm{max}$

We add the coverage via a COMSOL variable:

    comsol_variables['coverage']=['Kads*[[H2]]*Lmol/(1.+[[H2]]*Lmol*Kads)*theta_max','H Coverage according to Langmuir isotherm']

Note that the surface concentrations of species are indicated here by the double brackets [[...]]. Any species surface concentration can be referred like this.    

Finally, we can define the H$_2$ flux as:

    species['H2']['flux-equation'] = 
	    'rho*coverage*exp(-(Ga+alpha*F_const*(phiM-phi-phiEq))/RT)' #(mol/s/m^2)

We can also define a current density equation in the same way:

    species['H2']['current density-equation'] = 
	    '-2*F_const*rho*coverage*exp(-(Ga+alpha*F_const*(phiM-phi-phiEq))/RT)' #(mol/s/m^2)

Fixed flux expressions can be combined with flux-equation expressions and the remaining species fluxes will be automatically calculated. 

## CatMAP

The most advanced method of defining reactant fluxes is via a mean-field kinetic model. 
