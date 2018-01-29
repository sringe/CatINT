# How to define fluxes and reaction rates? 

## Fixed species fluxes 

Assume the HER in alkaline conditions:

    electrode_reactions={
        'H2': {'reaction': 'H2O + 2 e- -> H2 + 2 OH-'}
                        }
We can now define the flux of all the species e.g. by defining them in the species dictionary. For this particular reaction, this may look like this:

    species={
    	    'H2': { 	   'symbol': 'H_2',
		                   'name': 'hydrogen',
		                   'diffusion': 4.50e-009},
            'OH-': { 	   'symbol': 'OH^-',
		                   'name':                'hydroxide',
	                       'diffusion':            '5.273e-9',
	                       'bulk concentration':   1e-7},
	        }
We can add fluxes () by setting them individually:

    species['H2']['flux']=1e-4

Fluxes that are not given as input will be calculated from the other fluxes. In this case, the flux of OH<img src="https://rawgit.com/sringe/CatINT/master/svgs/bcc64d6c62f5ee4343ee672cf4fa895c.svg?invert_in_darkmode" align=middle width=10.274055pt height=26.17758pt/> can be calculated internally as:

<img src="https://rawgit.com/sringe/CatINT/master/svgs/03a16dd1a2fc2cd5c57a42dd8ca7ce95.svg?invert_in_darkmode" align=middle width=89.31549pt height=21.68331pt/>

So only one fluxe has to be defined per equation. 

## Rate-Equations

If multiple chemical reactions appear at the electrode, the flux of a species will be determined by many reactions at the same time 
