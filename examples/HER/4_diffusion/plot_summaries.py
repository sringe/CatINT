#!/usr/bin/env python
import pickle
import sys
import numpy as np
import os
from catmap.model import ReactionModel
import matplotlib.pyplot as plt
from matplotlib import rc
import matplotlib.mlab as mlab
from matplotlib.ticker import NullFormatter
from pylab import *
from matplotlib.mlab import griddata
from matplotlib.font_manager import FontProperties

class Object(object):
    pass

mathtext_prop = {'fontset' : 'custom',
                 'it' : 'serif:italic',
                 'sf' : 'Helvetica:bold',
                 'cal' : 'serif:italic:bold'}
rc('font', family='serif', serif='Helvetica')
rc('mathtext', **mathtext_prop)
rc('xtick', labelsize=15)
rc('ytick', labelsize=15)

def get_data(pickle_file):
    a = pickle.load(open(pickle_file))
    data = Object()
    #COVERAGES
    data.coverage_names = model.output_labels['coverage']
    coverage_map = np.array(a['coverage_map'])
    data.voltage = []    
    scaler_array = coverage_map[:,0]
    for s in scaler_array:
        data.voltage.append(s[0])
    coverage_mpf = coverage_map[:,1]
    data.coverage = np.zeros((len(coverage_mpf),len(data.coverage_names)))    
    for i in range(0,len(coverage_mpf)):
        for j in range(0,len(coverage_mpf[i])):
            float_rate = float(coverage_mpf[i][j])
            data.coverage[i][j]=float_rate
    #PRODUCT NAMES
    data.prod_names = model.output_labels['production_rate']
    production_rate_map = np.array(a['production_rate_map'])
    production_rate_mpf = production_rate_map[:,1]
    data.production_rate = np.zeros((len(production_rate_mpf),len(data.prod_names)))
    data.voltage = np.zeros((len(production_rate_mpf),1))
    for i in range(0,len(production_rate_mpf)):
        data.voltage[i][0] = production_rate_map[:,0][i][0]
        for j in range(0,len(data.prod_names)):
            float_rate = float(production_rate_mpf[i][j])
            data.production_rate[i][j]=float_rate
    #RATES
    data.rate_names = model.output_labels['rate']
    rate_map = np.array(a['rate_map'])
    rate_mpf = rate_map[:,1]
    data.rate = np.zeros((len(rate_mpf),len(data.rate_names)))
    for i in range(0,len(rate_mpf)):
        for j in range(0,len(rate_mpf[i])):
            float_rate = float(rate_mpf[i][j])
            data.rate[i][j]=float_rate
    return data

#SETTINGS
def convert_TOF(A): # Given a list, convert all the TOF to j(mA/cm2) using 0.161*TOF(According to Heine's ORR paper)
    B = [-0.161*rate for rate in A]
    return B
colormap = plt.get_cmap('gist_ncar')
color_list = ['black', 'darkgreen', 'blue', 'firebrick', 'darkviolet']

#PLOTS
pH = ['0','4','7','10','14']
for i in range(0,len(pH)): 
    j = pH[i]
    log_file = 'HER_pH'+j+'.log'
    model = ReactionModel(setup_file = log_file)
    pickle_file = 'HER_pH'+j+'.pkl'
    phX = get_data(pickle_file)
    idx=phX.prod_names.index('H2_g')
    data=np.column_stack((phX.voltage, phX.production_rate[:,idx]))
    plt.plot(data[np.argsort(data[:, 0])][:,0],convert_TOF(data[np.argsort(data[:, 0])][:,1]), color = color_list[i], linewidth=1.5, label='pH '+j)
    
fig1 = plt.figure(1, figsize=(7.5, 5.5))
ax = fig1.add_subplot(1,1,1)
leg = plt.legend(loc='lower right', ncol=1, prop={'size':12}, fancybox=True, shadow=False)
leg.get_frame().set_facecolor('none')
leg.get_frame().set_linewidth(0.0)
plt.gcf().subplots_adjust(bottom=0.18, left=0.18)
plt.xlim((-1.5,0.5))
plt.ylim((-12,1))
plt.xlabel(r'U (V vs. SHE)', fontsize=16)
plt.ylabel(r'i (mA/cm$^{2}$)', fontsize=16)
fig_name = 'polarization.png'
fig1.savefig(fig_name)
plt.close()

for i in range(0,len(pH)): 
    j = pH[i]
    log_file = 'HER_pH'+j+'.log'
    model = ReactionModel(setup_file = log_file)
    pickle_file = 'HER_pH'+j+'.pkl'
    phX = get_data(pickle_file)
    data=np.column_stack((phX.voltage, phX.coverage[:,0]))
    plt.plot(data[np.argsort(data[:, 0])][:,0],data[np.argsort(data[:, 0])][:,1], color = color_list[i], linewidth=1.5, label='pH '+j)

fig2 = plt.figure(1, figsize=(7.5, 5.5))
ax = fig2.add_subplot(1,1,1)
leg = plt.legend(loc='lower left', ncol=1, prop={'size':12}, fancybox=True, shadow=False)
leg.get_frame().set_facecolor('none')
leg.get_frame().set_linewidth(0.0)
plt.gcf().subplots_adjust(bottom=0.18, left=0.18)
plt.xlabel('coverage', fontsize=16)
plt.ylabel('U_SHE (V)', fontsize=16)
plt.xlim((-1.5,0.5))
#plt.ylim((-12,1))
plt.yscale('log')
plt.xlabel(r'U (V vs. SHE)', fontsize=16)
plt.ylabel(r'coverage', fontsize=16)
fig_name2 = 'coverages.png'
fig2.savefig(fig_name2)
plt.close()

for i in range(0,len(pH)): 
    j = pH[i]
    log_file = 'HER_pH'+j+'.log'
    model = ReactionModel(setup_file = log_file)
    pickle_file = 'HER_pH'+j+'.pkl'
    phX = get_data(pickle_file)
    data=np.column_stack((phX.voltage, phX.rate[:,2]))
    plt.plot(data[np.argsort(data[:, 0])][:,0],data[np.argsort(data[:, 0])][:,1], color = color_list[i], linewidth=1.5, label='pH '+j)

fig3 = plt.figure(1, figsize=(10.5, 5.5))
ax = fig2.add_subplot(1,1,1)
leg = plt.legend(loc='lower left', ncol=1, prop={'size':12}, fancybox=True, shadow=False)
leg.get_frame().set_facecolor('none')
leg.get_frame().set_linewidth(0.0)
plt.gcf().subplots_adjust(bottom=0.18, left=0.18)
plt.axvline(x=0.0, ymin=-10, ymax = 10, linewidth=2, color = color_list[0])
plt.axvline(x=-0.2355, ymin=-10, ymax = 10, linewidth=1, color = color_list[1])
plt.axvline(x=-0.4122, ymin=-10, ymax = 10, linewidth=1, color = color_list[2])
plt.axvline(x=-0.5888, ymin=-10, ymax = 10, linewidth=1, color = color_list[3])
plt.axvline(x=-0.824, ymin=-10, ymax = 10, linewidth=1, color = color_list[4])
plt.xlim((-1.2,0.1))
plt.ylim((-3.0,3.0))
plt.xlabel(r'U (V vs. SHE)', fontsize=16)
plt.ylabel(r'rate', fontsize=16)
plt.title('Tafel')
fig_name3 = 'rates.png'
fig3.savefig(fig_name3)
plt.close()


