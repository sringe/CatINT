===================================
CO2 reduction at polycrystalline Cu
===================================

Creation of the ``CatINT`` input file
-------------------------------------

``COMSOL`` arguments
~~~~~~~~~~~~~~~~~~~~

``COMSOL`` has to be run only once instead of iteratively as in the case of the ``CatMAP`` calculator. In that case, it is possible to let ``COMSOL`` do the iterations over the descriptor range, by setting

.. code:: python

    comsol_args['desc_method'] = 'internal'

Since the potential is correlated with the flux of species, it serves as a kind of non-linearity ramping and we do not need the ``'flux_factor'`` as an additional ramping. We thus define

.. code:: python

    comsol_args['par_name'] = 'phiM'
    comsol_args['par_method'] = 'internal-reinit'

Note, that ``'par_name'`` must be a name from the ``transport.descriptors`` list and the values are assigned automatically, if not changed here. We can perform the iterations by using the ``'internal-reinit'`` (solutions are reinitialized at each potential -- default) or ``'internal-cont'`` (solutions are initialized by the previous potential)  methods.
