tp_ref_data={\
#all constants & rates at room temperature
#Millero1997: http://www-naweb.iaea.org/napc/ih/documents/global_cycle/vol%20I/cht_i_09.pdf
#reactants:             (first line are the educts, second the products)
#constant:              (dimensionless (mol/m^3))
#rates:                 (forward and backward rates)
'electrolyte_reactions':
{
    'bicarbonate-base':
    {
        'buffer-base':      {   'reaction':            'CO2 + OH- <-> HCO3-',
                            #PURE WATER, Emerson
                            'constant':             43750.0,                                #Gupta:  44400.0 
                            'rates':                [7.0,16e-5]},                           #Gupta:  [5.93,13.4e-5]
                            #salinity S=35, Schulz2006
#                             'constant':             22966.014418, #Schulz2006, m^3/mol
#                             'rates':                [2.23,9.71e-5]}, #Schulz2006
    ##############################################################################################################
        'buffer-base2':     {   'reaction':            'HCO3- + OH- <-> CO32- + H2O', 
                            #PURE WATER ???? Gupta
                            'constant':              4.66,
                            'rates':                [1.0e5,21459.2274]},
                            #salinity S=35, Schulz2006
#                             'constant':             19.60784, #Schulz2006, m^3/mol
#                             'rates':                [6e6,306000]}, #Schulz2006
    },


    ##############################################################################################################
    ##############################################################################################################
    ##############################################################################################################


    'bicarbonate-acid':\
    {
    ##############################################################################################################
        'buffer-acid':      {   'reaction':            'CO2 + H2O <-> HCO3- + H+', 
                            #PURE WATER, Emerson
                            'constant':             0.000445,                               #Gupta:  0.000444 
                             'rates':               [3.7e-2,8.3333]},                       #
                            #salinity S=35, Schulz2006
#                             'constant':             0.00138951310,
#                             'rates':               [3.71e-2,26.7]},
    ##############################################################################################################
        'buffer-acid2':     {   'reaction':             'HCO3- <-> CO32- + H+',
                            #PURE WATER, Millero1997
                            #'constant':              4.79e-8,
                            'constant':             3.5317025629468759e-07,              #https://www.iaea.org/ocean-acidification/act7/Guide%20best%20practices%20low%20res.pdf
                            'rates':                [59.44,1.68304093e8]},                  #assuming Schulz2006 for hin-reactio !!!!!!!NOT SALINITY CORRECTED!!!!!!!
##                            'rates':                [59.44,12.409e8]},                  #assuming Schulz2006 for hin-reaction !!!!!!!NOT SALINITY CORRECTED!!!!!!!
                            #salinity S=35, Schulz2006
#                            'constant':             1.1888e-06,                            #Emerson: 1.0715e-6
#                            'rates':                [59.44,5e7]},                          
    ##############################################################################################################
    },

    ##############################################################################################################
    ##############################################################################################################
    ##############################################################################################################

    'phosphate-acid':\
        #from Ryu et al. Angew. Chemie (10.1002/anie.201802756)
        {
        'phosphate-1':  {   'reaction':     'H3PO4 + H+ <-> H2PO4-',
                            'constant':     0.00629*1000., #6.5e-3*1000,
                            'rates':        [5.6e8,8.9e10/1000.]
                        },
        'phosphate-2':  {   'reaction':     'H2PO4- + H+ <-> HPO42-',
                            'constant':     6.32e-8*1000., #6.2e-8*1000,
                            'rates':        [6.32e2, 1e10/1000.] #7.0,7/6.2e-8]
                        },
        'phosphate-3':  {   'reaction':     'HPO42- + H+ <-> PO43-',
                            'constant':     4.47e-13*1000., #3.6e-13*1000,                     #https://www.chem.fsu.edu/chemlab/Mastering/PhosphateBuffers.htm
                            'rates':        [4.47e-3, 1e10/1000.] #7.0,7/3.6e-13]                       #guessed
                        },
        },
    'citrate-acid':\
        #from Ryu et al. Angew. Chemie (10.1002/anie.201802756)
        #all estimated assuming diffusion limitations
        {
        'citrate-1':    {   'reaction':     'H3Cit + H+ <-> H2Cit-',
                            'constant':     0.000745*1000,
                            'rates':        [7.45e6,1e10/1000.]
                        },
        'citrate-2':    {   'reaction':     'H2Cit- + H+ <-> HCit2-',
                            'constant':     1.73e-5*1000,
                            'rates':        [1.73e5,1e10/1000.]
                        },
        'citrate-3':    {   'reaction':     'HCit2- + H+ <-> Cit3-',
                            'constant':     4.02e-7*1000,
                            'rates':        [4.02e3,1e10/1000.]
                        }
        },
    'borate-base':\
        {
        'borate-1':     {   'reaction':     'H3BO3 + OH- <-> H2BO3- + H2O',
                            'constant':     5.75e-10*1000,          #https://en.wikipedia.org/wiki/Boric_acid
                        },
        'borate-2':     {   'reaction':     'H2BO3- + OH- <-> HBO32- + H2O',
                            'constant':     3.98e-13*1000           #https://en.wikipedia.org/wiki/Boric_acid
                        },
        'borate-3':     {   'reaction':     'HBO32- + OH- <-> BO33- + H2O',
                            'constant':     5.01e-14*1000,          #https://en.wikipedia.org/wiki/Boric_acid
                        }
        },
                            
    ##############################################################################################################
    ##############################################################################################################
    ##############################################################################################################

    'water-diss':
        {'self-dissociation of water':            {   'reaction':             'H2O <-> OH- + H+',
                            'constant':             1e-8, #(mol/m^3)^2
                            'rates':                [2.4e-5*1000.,2.4e-5/1e-14/1000.]} #Singh
        }

    ##############################################################################################################
    ##############################################################################################################
    ##############################################################################################################
}
}
