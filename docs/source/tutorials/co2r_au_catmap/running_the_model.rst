Running the model
-----------------

Running the model, requires also to set-up a ``CatMAP`` mkm file and txt file containing all energies, see https://catmap.readthedocs.io/en/latest/ for details. Having created the ``CatINT`` input file ``run.py`` and the ``CatMAP`` files ``co2r.mkm`` and ``co2r.txt``, we can simply run ``CatINT`` via

.. code:: python

    python run.py

This will create a new folder preceding the name given in the definition of the ``transport`` instance and write all output files into this folder. 

``catmap_input``
~~~~~~~~~~~~~~~~

Contains subfolders ``desc_?`` for each descriptor value. This folder contains all input files and model definitions (also given in terms of free energy diagrams). 

``catmap_output``
~~~~~~~~~~~~~~~~~

Contains subfolders ``desc_?`` for each descriptor value. This folder contains all output files of ``CatMAP``, as coverages ``cov_[species].tsv``, reaction rates ``j_[species].tsv``, elementary step reaction rates ``jelem_?.tsv`` and degree of rate control ``rc_[product]_[species].tsv``

``comsol_results``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The folders ``comsol_results_id000_[index of descriptor 1]_[index of descriptor 2]`` contain all ``COMSOL`` results for a specific descriptor pair (In our case only a single descriptor, the potential, is used). The results outputted by ``COMSOL`` are defined in the ``transport.comsol_args['outputs']`` dictionary. The default outputs that are always generated are:

.. code:: python

    comsol_args['outputs']=['concentrations','electrostatics','electrode_flux','rho_charge']

These refer to species concentrations, electrostatics (potential and electric field), electrode fluxes and the charge density. Usually these outputs are enough, but if more are needed they can be easily asked for using the ``comsol_args['outputs']`` key in the ``CatINT`` input file and extending the ``catint/comsol_model.py`` if needed for the required output.

The folder also contains the compile and run log files for debugging purposes as well as copies of the ``COMSOL`` input file (``mph`` file) and ``.java`` input and compiled class ``.class`` file used for the particular calculations. Backups of the ``.java`` files with timestamps for debugging purpose are also saved in the results folder.

Convergence problems
~~~~~~~~~~~~~~~~~~~~

Convergence problems can happen in particular at high potentials. We found that the reason is often the very small CO2 concentration at the electrode which results in negative CO2 concentrations due to numerical issues. ``CatINT`` tries to solve this issue by disregarding the ``CatINT``-``COMSOL`` iteration step and using the previous CO2 concentration (but all updated species fluxes and surface concentrations). This mostly helps, but can also fail at higher potentials. Other possibilites of influencing convergence are a change of the ``grid_factor`` if ``COMSOL`` has problems converging.  A number of default steps are implemented in ``CatINT`` to try to get convergence, sometimes it even helps to restart the same ``COMSOL`` calculation due to slight numerical differences. If all does not help, ``CatINT`` stops and raises and error message.  If that is the case, we suggest to also try a finer potential descriptor range. The latter is important, since ``CatINT`` initializes the next descriptor step with the species surface concentrations of the previous step which yields a better initialization and avoid very bad first iteration steps which can break the convergence. Also the use of solver sequences might be possible to increase convergence. 
 

