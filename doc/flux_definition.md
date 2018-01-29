
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
    
The fluxes of the remaining species will be calculate automatically from the rate equation, i.e.

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

In COMSOL, we can define fluxes also as equations depending e.g. on the local concentrations of the species and the local potential. We can e.g. define the flux of H$_2$ in terms of a Butler-Volmer relation:

$j_\mathrm{H_2}=\rho\cdot\theta_\mathrm{H}\cdot A\cdot \exp(-[G_a+\alpha\cdot F\cdot(\Phi_\mathrm{M}-\Phi)]/RT)$,

where $\rho$ is the active site density, $\theta_\mathrm{H}$ is the coverage of H, $G_a$ is the activation barrier of the rate-determining step, $\alpha$ the transfer coefficient and $\Phi_\mathrm{M}-\Phi$ the difference between 

    species['H2']['flux-equation'] = 'rho*[[H2]]*exp(-Ga+(Vm-))' #(mol/s/m^2)
