.. _analysis:

Analysis tools
--------------

In general, ``CatINT`` includes two main tools for analyzing the results of the simulation:

- ``$CATINT/tools/plotting_catint.py``
- ``$CATINT/tools/plotting_catmap.py``

``plotting_catmap.py``
~~~~~~~~~~~~~~~~~~~~~~

``$CATINT/tools/plotting_catmap.py`` reads the results from the ``catmap_output`` folder and plots the current densities, coverages and degree of rate control. It supports also the comparison with experimental data, that is contained in the ``$CATINT/data/[reaction]/CSV`` folder.

``plotting_catint.py``
~~~~~~~~~~~~~~~~~~~~~~

``$CATINT/tools/plotting_catint.py`` reads the results from the ``comsol_results_?`` folder and plots all desired variables. Currently, ``CatINT`` supports a range of variables to be plotted. In general, for a variable named ``surface_[variable]``, the variable is plotted at the reaction plane  at a potential specified via ``--desc`` (maps to closest actually calculated potential). In all other cases, the variable is just plotted as a function of space (x):

- properties (``--prop``) supporting reaction plane (prepend ``surface_``) and spatial plot:
    - ``concentration``: Concentrations of all species
    - ``activity``: Activities of all species
    - ``activity_coefficient``: Activity coefficients of all species
    - ``efield``: Electric field
    - ``potential``: Electrostatic potential
    - ``charge_density``: Charge density
    - ``pH``: pH calculated using the proton activity (if in species list, otherwise the hydroxide activity)
- Other supported properties are:
    - ``electrode_current_density``: Partial current density of all products, should give the same results as the plotting of the ``catmap_output`` data using the ``$CATINT/tools/plotting_catmap.py`` script
    - ``pKw``: The water self-dissociation product at the reaction plane as a function of potential
    - ``concentration_at_x: Plot the concentrations of species at a specified distance x (specified using ``--xfixed``)

Also remember to set the potential scale correctly by specifying ``--scale RHE/SHE``.
