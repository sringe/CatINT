Creation of the ``CatINT`` input file
-------------------------------------

Import
~~~~~~

We start with the creation of the ``CatINT`` python input file which is contained in ``$CATINT/examples/02_CO2R_Au_CatMAP/run.py``. We import the main transport class and the calculator:

.. code:: python

    import numpy as np
    from catint.transport import Transport
    from catint.calculator import Calculator
    from catint.units import unit_NA

``CatINT`` comes with it's own units package.

System Settings
~~~~~~~~~~~~~~~

First, we define the system reaction conditions:

.. code:: python

    pH = 6.8
    system=\
        {
        #ENVIRONMENTAL CONDITIONS
        'temperature':  298,         #K
        'pressure':     1.013,       #bar
        'bulk_pH': pH,
        #MASS TRANSPORT
        'boundary thickness': 8.E-05, #m
        #ELECTROSTATICS
        'epsilon': 78.36,            #eps_0
        'migration': True,
        #REACTIONS
        'electrode reactions': True,
        'electrolyte reactions': True,
        #CHARGING
        'charging_scheme':'comsol',
        'phiM':-0.5,
        'phiPZC': 0.16,              #V vs. SHE
        'Stern capacitance': 20.,    #micro F/cm2
        #KINETICS
        'potential drop':'Stern',
        'active site density': 9.61-05/unit_NA*(1e10)**2, #mol sites/m^2
        #INITIALIZATION
       # 'init_folder':init_folder,
        }

``temperature`` and ``pressure`` are important for the finite temperature free energy corrections that are applied in ``CatMAP``, but the temperature also enters the mass transport equations. The ``pH`` should be adjusted to the reaction conditions. `boundary thickness`` refers to the boundary thickness which defines the cell dimensions. ``epsilon`` is the bulk dielectric constant that enters the mass transport equations. ``migration`` turns the migrational motion of species on or off. Further, electrode and electrolyte reactions can be turned on or off. The charging of the surface is controlled via the ``charging_scheme`` tag. Currently, ``CatINT`` supports different methods to define a surface charge density - potential relation :math:`\sigma(\phi^\mathrm{M})` which is used to evaluate potential-dependent free energies:

- ``'comsol'``:  Get the relation from ``COMSOL`` directly via the electrostatics defined there.
- ``'input'``:     Define a relation manually. This enables us to take the charging relation directly from the ``CatMAP`` input file. There, currently different ways are supported:
    - ``sigma_input=['CH',20.]``. Define surface charge density via :math:`\sigma=C_\mathrm{H} (\phi^\mathrm{M}-\phi^\mathrm{PZC})`, where :math:`C_\mathrm{H}` is a capacitance.
    - ``sigma_input='file.txt'``: Define surface charge density from  :math:`\sigma(\phi^\mathrm{M})` relation given in file name. The discrete data in the file is interpolated

For the electrostatics within the mass transport (PNP equations), the `Stern capacitance` (also called gap capacitance) is needed to define the Robin boundary condition. The same is valid for ``phiPZC``, :math:`\phi^\mathrm{M,PZC}` which is given on an SHE scale. ``phiM`` is the starting potential for the calculation which is usually overwritten when defining a descriptor range (see below). The `potential drop` defines how to calculate the driving electrochemical potential in ``CatMAP``, i.e. if a Frumkin correction is considered or not. If ``'potential_drop':'Stern'`` is selected a Frumkin correction is applied, so that the driving potential drop is :math:`\phi^\mathrm{M}-\phi^\ddagger` with :math:`\phi^\ddagger` is the potential at the reaction plane.

Finally ``CatINT`` supports to initialize a run form a previous one by using the ``init_folder`` tag. 

Electrolyte Reactions
~~~~~~~~~~~~~~~~~~~~~

We define electrolyte reactions which should be included in the mass transport. Currently ``CatINT`` supports multiple buffer reactions which are defined in the ``catint/data.py`` file. Here we include bicarbonate buffer reactions both using water and protons as donors (acidic and alkaline conditions) as well as the water self dissociation equilibria

.. code:: python

    electrolyte_reactions=['bicarbonate-base','water-diss',{'additional_cell_reactions':'bicarbonate-acid'}]

``CatINT`` evaluates concentrations at the boundary layer based on these buffer equilibria reactions. Since the acidic buffer reactions are already contained in the combination of alkaline and water dissociation, they are included as ``additional_cell_reactions``, meaning that they are active during the mass transport simulation, but not used for initializing the concentrations. 

Electrode Reactions
~~~~~~~~~~~~~~~~~~~

Now it is time to think about the reactions at the electrode, the electrochemical reactions. In our case, we consider the reduction of CO2 with water as a proton donor, which is defined as:

.. code:: python

    electrode_reactions={
        'CO': {'reaction': 'CO2 + H2O + 2 e- -> CO + 2 OH-'},
    }

Species Definitions
~~~~~~~~~~~~~~~~~~~

After defining the reactions, we need to also think about all species that our system, we define them as a dictionary:

.. code:: python

    species=\
        {
        'K+':               {'bulk_concentration':   'charge_neutrality',
                             'MPB_radius':           2*3.5e-10},
        'CO2':              {'bulk_concentration':   'Henry'},
        'OH-':              {'bulk_concentration':   10**(pH-14.)*1000.0}, #mol/m^3
        'CO':               {'bulk_concentration':   0.0}
        }
    
All species that are part of the electrolyte reactions are already automatically added to the species dictionary in ``CatINT``. All other species have to be added here. Charges are automatically assigned according to the species name. We add potassium cations which should neutralize all anions so that the system is charge neutral in the bulk solution (boundary layer). We chose a ``MPB_radius`` of 3.5 Angstrom which is important since the negative electrode potental dramatically increases the potassium concentrations. We then define the CO2 concentration at the boundary layer to be given by Henry's law which will use the pressure defined in system settings and the Henry constant in ``data/henry_constants.txt`` to evaluate the equilibrium CO2 concentration. The buffer component concentrations are now evaluated using the equilibrium buffer equations and must not be specified. It is also possible though to specify a buffer concentration and let ``CatINT`` calculate the CO2 concentration via the buffer equations. Finally, we set the concentration of hydroxide anions according to the pH (proton concentration are automatically evaluated using the water dissociation equilibrium). All concentrations are given in :math:`\mathrm{mol}/\mathrm{m}^3`.

Descriptors
~~~~~~~~~~~

In a common application, ``CatINT`` calculations should be run for a specified parameter or descriptor range. In this example, we want to simulate a polarization curve, our descriptor is therefore the electrode potential :math:`\phi^\mathrm{M}`, we define it like this:

.. code:: python

    phimin=-0.5
    phimax=-2.0
    dphi=0.01
    descriptors={'phiM':list(np.linspace(phimin,phimax,-(phimax-phimin)/dphi+1))}

Currently ``CatINT`` supports only the potential as descriptor, others could be implemented, if needed. The ``CatINT`` calculator iterates over the descriptor list and solves the coupled mass transport -- microkinetics model at each potential.

``COMSOL`` arguments
~~~~~~~~~~~~~~~~~~~~

There are a couple of predefined ``COMSOL`` arguments which are saved in the ``tp.comsol_args`` dictionary (suppose that ``tp`` is the transport instance that we create in the end of this tutorial page). We can define COMSOL variables and parse them to CatINT via the ``comsol_args`` tag:

.. code:: python

    comsol_args={}
    #parameter
    comsol_args['parameter']={}   
    comsol_args['parameter']['grid_factor']=[str(100),'Grid factor']
    comsol_args['parameter']['grid_factor_domain']=[str(100),'Grid factor for domain']
    comsol_args['parameter']['grid_factor_bound']=[str(200),'Grid factor for boundary']
    #solver_settings
    comsol_args['solver_settings']={}
    comsol_args['solver_settings']['direct']={}
    comsol_args['solver_settings']['direct']['nliniterrefine']=True
    comsol_args['solver_settings']['ramp']={}
    comsol_args['solver_settings']['ramp']['names']=['PZC','CS']
    comsol_args['solver_settings']['ramp']['dramp']=0.01
    #par_method
    comsol_args['par_method']='internal'

    #SOLVER SEQUENCE
    #comsol_args['solver_settings']['solver_sequence']='tds_elstat'
    #OTHER PARAMETER
    #comsol_args['parameter']['RF']=[1,'Roughness Factor']

This is in particular useful for modifying numerical solver settings. In our case, we first define a ``grid_factor`` which tells ``COMSOL`` about the minimal finite element mesh width. A higher factor means a finer mesh and the mesh can be defined for the domain and boundary separately. Parameter definitions always a list of two entries, the value parsed as a string and the name or description inside COMSOL. Some more specific numerical parameters can be edited and changed here to help convergence. In particular, the ``'ramp'`` flag enables to slowly ramp non-linearities in the equations, in our case it slowly ramps up the PZC and the Helmholtz/Stern/gap capacitance which can be useful if the systems has a PZC far from the initial potential to be evaluated. A flux ramping is always applied and controlled by the `'dramp'` flag which defines the interval in which the fluxes are ramped from 0 to 100 %.  Additional possible settings involve the definition of solver sequences to improve convergence (e.g. first solving electrostatics only, then the coupled mass transport/electrostatic problem). Also, it is possible to define a roughness factor that multiplies all fluxes by a constant.

Inside ``CatINT``, a couple of ``COMSOL`` variables are assigned by default. The iterations over ``tp.descriptors`` are performed inside ``CatINT`` and ``COMSOL`` is compiled and run for each descriptor value. This behavior is defined via the ``COMSOL`` key:

.. code:: python

    comsol_args['desc_method'] = 'external'

Inside ``COMSOL``, we have the possibility to sweep over a particular parameter space to enable convergence. A common way to do this, is to ramp up the non-linearities in the equations as the flux of the species. This is the default in ``CatINT``:

.. code:: python

    comsol_args['par_name'] = 'flux_factor'
    comsol_args['par_values'] = 'range(0,'+str(comsol_args['solver_settings']['ramp']['dramp'])+',1)'
    comsol_args['par_method'] = 'internal'

The range can be modified by the ``dramp`` key as discussed before. The ``par_method`` key indicates the way that ``COMSOL`` should treat the parametric sweep: ``'internal'`` uses an auxiliary parameter sweep, while ``'external'`` uses a regular parameter sweep. Although both should do in principle the same, there fine differences between both methods inside ``COMSOL``, and generally the ``'internal'`` sweep seems to be more stable.

``CatMAP`` arguments
~~~~~~~~~~~~~~~~~~~~

Some additional arguments can be parsed to the ``CatMAP`` calculator:

.. code:: python

    catmap_args={}
    #CATMAP DESCRIPTOR RAMPING
    catmap_args['desc_method']='automatic'
    #catmap_args['min_desc']=0.0
    catmap_args['min_desc_delta']=0.2
    catmap_args['max_desc_delta']=0.2
    #INTERACTIONS
    catmap_args['n_inter']='automatic'

In a regular ``CatMAP``-``COMSOL`` iteration loop, a single ``CatMAP`` calculation is required at the descriptor value of choice. This is referred to as ``desc_method='single_point'``. Sometimes, however, some potential values are hard to converge and it is better to provide ``CatMAP`` with a range of potentials. ``CatMAP`` will then try to find stable starting points and solve for all potentials and ``CatINT`` selects the potential that was actually needed. This happens if we choose ``desc_method='automatic'``. For this case, ``min_desc`` specifies the minimum descriptor value in the new created list of descriptors. Alternatively, ``min/max_desc_delta`` can be used to create a new list of descriptors around the current descriptor value. The descriptor range can be also taken from the ``CatMAP`` input file (``desc_method='from_input'``).

Flux definition
~~~~~~~~~~~~~~~

Now it is time to define how fluxes are evaluated within the ``CatINT`` model. In our case, we will use ``CatMAP`` to define fluxes but multiple options are available (cf. :ref:`flux-definition`).

.. code:: python

    species['CO']['flux']='catmap' #CO production rate
    species['CO2']['flux']='catmap' #CO2 consumption rate

Defining transport instance and assigning calculator
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

We are now ready to define the transport instance to which we parse all the previously defined dictionaries:

.. code:: python

    nx=200
    tp=Transport(
        species=species,
        electrode_reactions=electrode_reactions,
        electrolyte_reactions=electrolyte_reactions,
        system=system,
        catmap_args=catmap_args,
        comsol_args=comsol_args,
        model_name='CO2R',
        descriptors=descriptors,
        nx=nx)

``nx`` is hereby the starting number of finite elements. By using the ``COMSOL`` calculator this is usually relevant, because the ``'grid_factor_?'`` keys will define the mesh. However, ``tp.nx`` will be updated and thus the final size of the mesh within ``COMSOL`` can be accessed via this.

We now choose a calculator for the transport instance, in our case the ``COMSOL`` calculator and then assign the transport instance to a calculator object.

.. code:: python

   tp.set_calculator('comsol')
   c=Calculator(transport=tp,tau_scf=0.008,ntout=1,dt=1e-1,tmax=10,mix_scf=0.02)

Now, we can run the calculation:

.. code:: python

   c.run()

Mass transport-free CatMAP simulation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

We can also use ``CatINT`` just to run ``CatMAP``. This can be useful for easily comparing non-mass transport corrected and mass transport corrected results, or to just use the analysis tools of ``CatINT`` for ``CatMAP`` (cf. :ref:analysis). In order to do this, instead of assigning the ``tp`` instance to a ``Calculator`` object, we define a ``CatMAP`` instance and then run ``CatMAP`` for each potential:

.. code:: python

   cm=CatMAP(transport=tp,model_name='CO2R')
   for pot in descriptors['phiM']:
      tp.descriptors['phiM']=[pot]
      cm.run()

Mass transport extrapolation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Sometimes, ``COMSOL`` does not converge any more, but we have sufficient data that we think we can try to extrapolate all species concentrations at the reaction plane to a different potential region. This is done by polynomial functions, but these can be in principle specified by the user. We first import the ``extrapolate`` function:

.. code:: python

   from tools.extrapolate_surface_conc import extrapolate

We then start by running a mass-transport corrected ``CatINT`` simulations for the region of potentials which converges. Then we run the input file again until the assignment of a calculator:

.. code:: python

   tp.set_calculator('comsol')
   extra=extrapolate(tp=tp,extrapol_folder='CO2R_results_to_extrapolate')
   extra.plot()

This first initializes the ``tp`` instance as before, but that uses the old ``CatINT`` results folder which we named ``CO2R_results_to_extrapolate`` here to extrapolate concentrations of species at the reaction plane. The plot function enables to visualize the extrapolation and play around with the extrapolation functions.

If one is satisfied, one can remove the two last lines and instead define a ``CatMAP`` only calculator and use the extrapolated surface concentrations at each potential for the ``CatMAP`` calculation:

.. code:: python

   tp.set_calculator('comsol')
   cm=CatMAP(transport=tp,model_name='CO2R')
   for pot in descriptors['phiM']:
      for sp in tp.species:
         tp.species[sp]['surface_concentration']=10**extra.extrapol_func[sp](pot)
      tp.system['potential']=[extra.extrapol_func['voltage_diff_drop'](pot)]
      tp.system['surface_pH']=extra.extrapol_func['surface_pH'](pot)
   cm.run()

It is important to also extrapolate the ``voltage_diff_drop``, the Frumkin correction for the driving force as well as the pH at the reaction plane which both enter the ``CatMAP`` kinetics.
