#!/usr/bin/env python


from numpy import arange,array,ones,linalg
from matplotlib.pyplot import plot,show
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import rc
import re
from itertools import cycle
from units import *

import sys
import csv
import sys

class Object(object):
    pass

font = {'family' : 'serif',
        'serif':['Times'],
         'size'   : 16}
rc('font', **font)
rc('text', usetex=True)   
rc('xtick', labelsize=16)
rc('ytick', labelsize=16)

X=7
figsize = (3*X,2*X)

color = ['r','darkorange','gold','g','c','b','violet','r','darkorange','gold','g','c','b','violet']*2
linestyle = [':',':','-','-','-','-','-','-','-','-','-','-','-','-','-']
marker=['o','o','o','o','o','o','o','v','v','v','v','v','v','v']



def get_data(filename):
    o = open('CSV/'+filename,'rU')
    DATA = Object()
    DATA.data = np.array(list(csv.reader(o)))
    DATA.label=[str(d) for d in DATA.data[0,:]] #data[0,1:4]
    DATA.data = np.delete(DATA.data,0,0) # delete header
    return DATA


DATA_hs=get_data('COR_hori_ph_normalized_SHE.csv')
DATA_hr=get_data('COR_hori_ph_normalized_RHE.csv') #hori_data_RHE.csv')
#DATA_js=get_data('jaramillo_kanan_data_she.csv')
#DATA_jr=get_data('jaramillo_kanan_data_rhe.csv')

DATA_jr=get_data('CO2R_jaramillo_normalized_rhe.csv')

#DATA_jr2=get_data('CO2R_jaramillo_normalized_rhe.csv')

#DATA_h2s=get_data('hori_jpart_CO2R_SHE.csv')
DATA_h2s=get_data('CO2R_hori_normalized_SHE.csv')

DATA_kr=get_data('CO2R_kanan_normalized_RHE.csv')

#DATA_ss=get_data('strasser_jpart_CO2R_SHE.csv')
DATA_ss_005=get_data('strasser_jpart_from_prate_KHCO3_005_CO2R_SHE.csv')
DATA_ss_02=get_data('strasser_jpart_from_prate_KHCO3_02_CO2R_SHE.csv')
DATA_ss_01=get_data('strasser_jpart_from_prate_KHCO3_01_CO2R_SHE.csv')
DATA_ss_005_FE=get_data('strasser_jpart_KHCO3_005_CO2R_NHE.csv')
DATA_ss_01_FE=get_data('strasser_jpart_KHCO3_01_CO2R_NHE.csv')
DATA_ss_02_FE=get_data('strasser_jpart_KHCO3_02_CO2R_NHE.csv')

DATA_lr=get_data('COR_lei_normalized_rhe.csv')

jaramillo_mass_transport_limit = [4,4,4,4,2,2,0,0,4,4,2,2,2,2]

# o = open('hori_data_SHE.csv','rU')
# data = np.array(list(csv.reader(o)))
# DATA = Object()
# DATA.label=[str(d) for d in data[0,:]] #data[0,1:4]
# data = np.delete(data,0,0) # delete header


#for j in range(0,len(DATA.label)):
#    if j % 2 == 1:

#def tafel_fitting(x,y):
#    #fit a 2nd order polynomial first
#    ref=np.polyfit(x,y,2)
#    refd=np.poly1d(ref)
#    
#    #now iteratively remove single points and check when the slope is not changing much anymore
#    i=-1
#    for xx,yy in zip(x,y):
#        i+=1
#        fit=np.polyfit(x[i:],y[i:],1)
#        #compare
#        ref[0]
#        fitd=np.poly1d(
#

def plot_stuff(list_of_data,DATA,fit,joinlines, skip={}, voltage_mode='previous',take_log=False,linestyle=None,symbol=None,convert=None,fit_tafel=False):
    #if voltage_mode='previous' take previous point before list_of_data as voltages, if 'first', take first column for all data
    #if fit=0 Tafel, fit=1 poly, fit=2 poly no points
    #fit_tafel: New tafel fitting based on error estimation due to curvature
    n=0
    global symbols
    for j in list_of_data: #[1,3,5,7,9,11,13]:
        if symbol is not None and j==0:
            cmarker=symbol
        elif symbol is not None:
            cmarker=next(symbols)
        voltage=[]
        current=[]
        for k in range(0,len(DATA.data[:,j-1])):
            try:
                current.append(float(DATA.data[k,j]))
                if voltage_mode=='previous':
                    voltage.append(float(DATA.data[k,j-1]))
                elif voltage_mode=='first':
                    voltage.append(float(DATA.data[k,0]))
                else:
                    print 'Invalid voltage_mode'
            except ValueError:
                pass
        for k in range(0,len(voltage)):
            voltage_all.append(voltage[k])
            if take_log:
                current[k]=np.log10(current[k])
            current_all.append(current[k])
        linestyle_2 = "None"
        if joinlines:
            linestyle_2=":"
        if linestyle is None:
            linestyle = "None"
            if joinlines:
                linestyle=":"
        if any([a in DATA.label[j] for a in ['n-PrOH']]):
            color[n]='lightblue'
        elif any([a in DATA.label[j] for a in ['C$_2$H$_4$','CH$_3$COO','EtOH','C$_{2+}$','C$_2$']]):
            color[n]='r'
        elif any([a in DATA.label[j] for a in ['CH$_4$','C$_1$']]):
            color[n]='orange'
        elif 'HCOO' in DATA.label[j]:
            color[n]='b'
        elif 'H$_2' in DATA.label[j]:
            color[n]='k'
        elif 'CO' in DATA.label[j]:
            color[n]='olive'
        else:
            pass
            #color[n]='dark red'
        val=re.findall('pH[ ]*=[ ]*(\d+)',DATA.label[j])
        if len(val)>0:
            pH=float(val[0])
        else:
            pH=7.0
        if symbol is None:
            cmarker=marker[n]
        if convert is not None:
            if convert.split('_')[0]=='RHE':
                shift=-0.059*pH
            else:
                shift=+0.059*pH
        else:
            shift=0.0
        voltage=[v+shift for v in voltage]

        #label
        #labels.append(DATA.label[j])
        #legend_labels = numpy.concatenate([label_row_1, label_j_1, label_empty * 3, label_j_2, label_empty * 3, label_j_3, label_empty * 3])

        #if fit_tafel:
        #    tafel_fitting(voltage,current)
        print 'FIT',fit
        if not fit>2:
            plt.plot(voltage, current, color=color[n], linestyle=linestyle_2, marker = cmarker, label=DATA.label[j])  #linestyle = ':',
        if fit>0:
            if len(voltage)>1:
                species=DATA.label[j].split(',')[0].strip()
                skip_val=0
                for sk in skip:
                    if sk in species:
                        skip_val=skip[sk]
                if fit==1:
                    plot_tafel(voltage,current,color[n],skip_val,linestyle=linestyle)
                if fit>1:
                    V=voltage
                    J=current
                    linecol=color[n]
                    #if skip:
                    #    V=V[skip:]
                    #    J=J[skip:]
                    if 'H$_2' in DATA.label[j]:
                        order=3
                    else:
                        order=2
                    z=np.polyfit(V,J,order)
                    p=np.poly1d(z)
                    xfine=np.linspace(min(V),max(V),1000)
                    if linestyle=="None":
                        linestyle='-'
                    if fit>2:
                        plt.plot(xfine,p(xfine),linestyle=linestyle,color=linecol,label=DATA.label[j])
                    else:
                        print 'plotting fit',DATA.label[j],linestyle
                        print 'data',zip(V,J,p(V))
                        plt.plot(xfine,p(xfine),linestyle=linestyle,color=linecol)
                        
        n=n+1   
    plt.ylabel('log$j$ [mA/cm$^2_{\mathrm{real}}$] ')
    plt.ylim((-3,1))
    plt.xlim((-1.4,-1.05))
    leg = plt.legend(loc=1, 
          ncol=1, fontsize=10, numpoints=1) #, fancybox=True, shadow=False, loc='upper center', bbox_to_anchor=(0.5, 1.05))
    for line,text in zip(leg.get_lines(), leg.get_texts()):
        text.set_color(line.get_color())


def plot_eqm():
    plt.axvline(x=0.17, color='g')
    plt.text(0.185,1,'$E^o$(CH$_4$)', rotation='vertical', color='g', fontsize=13)
    plt.axvline(x=0.06, color='m')
    plt.text(0.02,1,'$E^o$(C$_2$H$_4$)', rotation='vertical', color='m', fontsize=13)
    plt.axvline(x=0.08, color='k')
    plt.text(0.095,1,'$E^o$(C$_2$H$_5$OH)', rotation='vertical', color='k', fontsize=13)

def plot_tafel(V,J,linecol,skip=0,linestyle='-'):
    # fit lines
    data=[V,J]
    data=map(list,zip(*data))
    data=np.array(data)
    print 'before',data
    #print 'before',data
    #data=sorted(data,key=lambda l:l[1], reverse=True)
    #data=np.array(data)
    #print 'after',data
    data.view(data.dtype.str+','+data.dtype.str).sort(order=['f0'], axis=0)
    print 'after',data
    V=data[:,0]
    J=data[:,1]
    if skip:
        V=V[skip:]
        J=J[skip:]
    A = array([V, ones(len(V))])
    w = linalg.lstsq(A.T,J )[0]
    line_tafel = w[0]*np.array(V)+w[1]
    tafel_slope = 1000/w[0]
    plt.plot(V,line_tafel,color=linecol, linestyle=linestyle)
    plt.text(V[0],line_tafel[0],str(int(tafel_slope))+'mV/dec', color=linecol, fontsize=13)  #

def plot_data(reference=['all'],species=['all'],pH=['all'],ci_bic=['all'],scale='RHE',reaction='all',system=['all'],coloring='species',fit_tafel=False):
    """
    ----------------------------------------------------------------------
    Wrapper to plot any data set of interest
    ----------------------------------------------------------------------
    reference: all, jaramillo, lei, kanan, strasser or hori (default: all)
    species:   all, C1-sum, C1, C2-sum, C2, H$_2$, CH$_4$, C$_2$H$_4$, ... [a list!]
    pH:        list of pH's to plot [list]
    ci_bic:    all, 0.05, 0.1, 0.2 [list of concentraitions of bicarbonate]
    scale:     RHE or SHE
    reaction:  all, CO2R or COR [list]
    system:    all, pc-Cu, OD-Cu, Cu [list]
    coloring:  which descriptor should be used to color results: species, pH, 
    fit_tafel: fit a tafel line to the low overpotential points until curvature is too strong (transport limitations)
    """

    global voltage_all
    global current_all
    global filter_reaction
    global filter_system

    filter_reaction=reaction
    filter_system=system

    voltage_all=[]
    current_all=[]

    if any([type(a) != list for a in [reference,species,pH,ci_bic,system]]):
        print 'Reference, species and pH must be lists!'
        sys.exit()

    #work on the species list, expand C1 and C2 in the respective species:
    #C1 does by definition here NOT include HCOO (formate,formic acid)!!
    c1_list=['CH$_4$','MeOH','CH$_3$OH'] #,'HCOO','MeOH','CH$_3$OH'] #,'C$_1$']
    c2_list=['C$_2$H$_4$','CH$_3$COO','EtOH','CH$_2$CH$_2$','C$_2$H$_6$','MeCHO','GlycAld','AcetAld','EtGlyci']
    c2p_list=c2_list+['n-PrOH','PrOH','AllylAlc','Acetone']
    species_reduced=[]
    for sp in species:
        if not ('C1' in species and sp in c1_list) and not ('C2' in species and sp in c2_list) and not ('C2+' in species and sp in c2p_list)\
                and not sp=='C1' and not sp=='C2' and not sp=='C2+':
            species_reduced.append(sp)
    if 'C1' in species:
        species_reduced+=c1_list
    if 'C2' in species:
        species_reduced+=c2_list
    if 'C2+' in species:
        species_reduced+=c2p_list
    species=species_reduced

    species=['C$_1$' if sp=='C1-sum' else sp for sp in species]
    species=['C$_2$,' if sp=='C2-sum' else sp for sp in species]
    species=['C$_{2+}$' if sp=='C2+-sum' else sp for sp in species]

    print('Plotting species = {}'.format(species))

    plt.figure()

    def isref(ref):
        if 'all' in reference:
            return True
        elif ref in reference:
            return True
        else:
            return False

    def find(search,large):
        f=re.findall(search,large)
        if len(f)>0:
            return f[0]
        else:
            return None

    def s2i(species,DATA,pH='7.0'):

        def anyph(pH,line):
            if pH=='all':
                return True
            found=re.findall('pH = (\d*\.?\d+)',line)
            if len(found)>0:
                pH_label=float(found[0])
            else:
                print 'No pH given in data file ',DATA.label
                print 'This is the line that was searched:',line
                sys.exit()
            if abs(pH_label-float(pH))<1e-3:
                return True
            else:
                return False

        def reac(line):
            if filter_reaction=='all':
                return True
            found=re.findall('('+filter_reaction+')',line)
            if len(found)>0:
                return True
            else:
                return False

        def system(line):
            for system in filter_system:
                if system=='all':
                    return True
                print 'checking',system,line
                found=re.findall(system,line)
                if len(found)>0:
                    return True
            return False

        heads=None
        if 'all' in species:
            heads=[]
            for i,h in enumerate(DATA.label):
                if i!=0:
                    if len(h)>0 and anyph(pH,h) and reac(h) and system(h) and i not in heads:
                        heads+=[i]
        else:
            heads=[]
            for sp in species:
                for i,h in enumerate(DATA.label):
                    print sp, h
                    if h.startswith(sp) and anyph(pH,h) and reac(h) and system(h) and i not in heads:
                        heads+=[i]
        return heads

    linestyles=cycle(['-','-',':','-.'])
    global symbols
    symbols=cycle(['d','D','o','x','p','*','v','h','1','2'])

    if fit_tafel:
        fit=1
    else:
        fit=2
    
    print DATA_ss_01.label

    data_labels=[','.join(a.label[1].split(',')[1:]) for a in [DATA_ss_01,DATA_ss_005,DATA_ss_02,DATA_ss_005_FE,DATA_ss_01_FE,DATA_ss_02_FE,DATA_jr,DATA_kr,DATA_lr,DATA_h2s,DATA_hr,DATA_hs]]

    print data_labels



    for cpH in pH:
        skip_dict={}
        for label in data_labels:
            if 'Hori' in label:
                skip_dict[label]={'HCOO':   5,\
                                  'CO':     9,\
                                  'H$_2$':  10}
            elif 'Jaramillo' in label:
                skip_dict[label]={'HCOO':   6,\
                                  'CH$_4$': 2,\
                                  'H$_2$':  5,\
                                  'C$_2$H$_4$':3,\
                                  'EtOH':3}
            elif 'Kanan' in label:
                skip_dict[label]={'HCOO':   9,\
                                  'CO':     10,\
                                  'H$_2$':     9,\
                                  'C$_2$H$_4$':2,\
                                  'C$_2$H$_6$':3}
            elif 'Strasser' in label:
                skip_dict[label]={'HCOO':   0,\
                                  'CO':     8}
        if abs(float(cpH)-13)<1e-5:
            for label in data_labels:
                if 'Jaramillo' in label:
                    skip_dict[label]['CH$_4$']=0
        for bic in ci_bic:
            if isref('strasser'):
                if bic in ['0.1','all']:
                    DATA=DATA_ss_01
                    name=','.join(DATA.label[1].split(',')[1:])
                    print 'name',name
                    skip=skip_dict[name]
                    spp=s2i(species,DATA,pH=cpH)
                    if len(spp)>0:
                        linestyle=next(linestyles)
                        symbol=next(symbols)
                        if scale=='RHE':
                            plot_stuff(spp,DATA,fit,0,skip,linestyle=linestyle,symbol=symbol,take_log=True,convert='SHE_TO_RHE',fit_tafel=fit_tafel)
                        elif scale=='SHE':
                            plot_stuff(spp,DATA,fit,0,skip,linestyle=linestyle,symbol=symbol,take_log=True,fit_tafel=fit_tafel)
                    #alternatively from FE's
                    DATA=DATA_ss_01_FE
                    name=','.join(DATA.label[1].split(',')[1:])
                    skip=skip_dict[name]
                    spp=s2i(species,DATA,pH=cpH)
                    if len(spp)>0:
                        linestyle=next(linestyles)
                        symbol=next(symbols)
                        if scale=='RHE':
                            plot_stuff(spp,DATA,fit,0,skip,linestyle=linestyle,symbol=symbol,voltage_mode='first',take_log=True,convert='SHE_TO_RHE',fit_tafel=fit_tafel)
                        else:
                            plot_stuff(spp,DATA,fit,0,skip,linestyle=linestyle,symbol=symbol,voltage_mode='first',take_log=True,fit_tafel=fit_tafel)
                if bic in ['0.05','all']:
                    DATA=DATA_ss_005
                    name=','.join(DATA.label[1].split(',')[1:])
                    skip=skip_dict[name]
                    spp=s2i(species,DATA,pH=cpH)
                    if len(spp)>0:
                        linestyle=next(linestyles)
                        symbol=next(symbols)
                        if scale=='RHE':
                            plot_stuff(spp,DATA,fit,0,skip,linestyle=linestyle,symbol=symbol,take_log=True,convert='SHE_TO_RHE',fit_tafel=fit_tafel)
                        elif scale=='SHE':
                            plot_stuff(spp,DATA,fit,0,skip,linestyle=linestyle,symbol=symbol,take_log=True,fit_tafel=fit_tafel)
                    #alternatively from FE's
                    DATA=DATA_ss_005_FE
                    name=','.join(DATA.label[1].split(',')[1:])
                    skip=skip_dict[name]
                    spp=s2i(species,DATA,pH=cpH)
                    if len(spp)>0:
                        linestyle=next(linestyles)
                        symbol=next(symbols)
                        if scale=='RHE':
                            plot_stuff(spp,DATA,fit,0,skip,linestyle=linestyle,symbol=symbol,take_log=True,voltage_mode='first',convert='SHE_TO_RHE',fit_tafel=fit_tafel)
                        else:
                            plot_stuff(spp,DATA,fit,0,skip,linestyle=linestyle,symbol=symbol,voltage_mode='first',take_log=True,fit_tafel=fit_tafel)
                if bic in ['0.2','all']:
                    DATA=DATA_ss_02
                    name=','.join(DATA.label[1].split(',')[1:])
                    skip=skip_dict[name]
                    spp=s2i(species,DATA,pH=cpH)
                    if len(spp)>0:
                        print 'taking indices',spp
                       # sys.exit()
                        linestyle=next(linestyles)
                        symbol=next(symbols)
                        if scale=='RHE':
                            plot_stuff(spp,DATA,fit,0,skip,linestyle=linestyle,symbol=symbol,take_log=True,convert='SHE_TO_RHE',fit_tafel=fit_tafel)
                        elif scale=='SHE':
                            plot_stuff(spp,DATA,fit,0,skip,linestyle=linestyle,symbol=symbol,take_log=True,fit_tafel=fit_tafel)
                    #alternatively from FE's
                    DATA=DATA_ss_02_FE
                    name=','.join(DATA.label[1].split(',')[1:])
                    skip=skip_dict[name]
                    spp=s2i(species,DATA,pH=cpH)
                    if len(spp)>0:
                        linestyle=next(linestyles)
                        symbol=next(symbols)
                        if scale=='RHE':
                            plot_stuff(spp,DATA,fit,0,skip,linestyle=linestyle,symbol=symbol,take_log=True,voltage_mode='first',convert='SHE_TO_RHE',fit_tafel=fit_tafel)
                        else:
                            plot_stuff(spp,DATA,fit,0,skip,linestyle=linestyle,symbol=symbol,take_log=True,voltage_mode='first',fit_tafel=fit_tafel)
            if bic in ['0.1','all']:
                if isref('jaramillo'):
                    DATA=DATA_jr
                    name=','.join(DATA.label[1].split(',')[1:])
                    skip=skip_dict[name]
                    spp=s2i(species,DATA,pH=cpH)
                    if len(spp)>0:
                        linestyle=next(linestyles)
                        symbol=next(symbols)
                        if scale=='RHE':
                            plot_stuff(spp,DATA,fit,0,skip,linestyle=linestyle,symbol=symbol,voltage_mode='first',take_log=True,fit_tafel=fit_tafel)
                        elif scale=='SHE':
                            plot_stuff(spp,DATA,fit,0,skip,linestyle=linestyle,symbol=symbol,voltage_mode='first',take_log=True,convert='RHE_TO_SHE',fit_tafel=fit_tafel)
                    DATA=DATA_lr
                    name=','.join(DATA.label[1].split(',')[1:])
                    skip=skip_dict[name]
                    spp=s2i(species,DATA,pH=cpH)
                    if len(spp)>0:
                        linestyle=next(linestyles)
                        symbol=next(symbols)
                        if scale=='RHE':
                            plot_stuff(spp,DATA,fit,0,skip,linestyle=linestyle,symbol=symbol,voltage_mode='first',take_log=True,fit_tafel=fit_tafel)
                        else:
                            plot_stuff(spp,DATA,fit,0,skip,linestyle=linestyle,symbol=symbol,voltage_mode='first',take_log=True,convert='RHE_TO_SHE',fit_tafel=fit_tafel)
                if isref('hori'):
                    DATA=DATA_h2s
                    name=','.join(DATA.label[1].split(',')[1:])
                    skip=skip_dict[name]
                    spp=s2i(species,DATA,pH=cpH)
                    if len(spp)>0:
                        linestyle=next(linestyles)
                        symbol=next(symbols)
                        if scale=='RHE':
                            plot_stuff(spp,DATA,fit,0,skip,voltage_mode='first',take_log=True,linestyle=linestyle,symbol=symbol,convert='SHE_TO_RHE',fit_tafel=fit_tafel) #,skip=4) #
                        elif scale=='SHE':
                            plot_stuff(spp,DATA,fit,0,skip,voltage_mode='first',take_log=True,linestyle=linestyle,symbol=symbol,fit_tafel=fit_tafel) #,skip=4) #

                    if scale=='RHE':
                        DATA=DATA_hr
                        name=','.join(DATA.label[1].split(',')[1:])
                        skip=skip_dict[name]
                    elif scale=='SHE':
                        DATA=DATA_hs
                        name=','.join(DATA.label[1].split(',')[1:])
                        skip=skip_dict[name]

                    spp=s2i(species,DATA,pH=cpH)
                    if len(spp)>0:
                        linestyle=next(linestyles)
                        symbol=next(symbols)
                        plot_stuff(spp,DATA,fit,0,skip,linestyle=linestyle,symbol=symbol,fit_tafel=fit_tafel) 

                if isref('kanan'):
                    DATA=DATA_kr
                    name=','.join(DATA.label[1].split(',')[1:])
                    skip=skip_dict[name]
                    spp=s2i(species,DATA,pH=cpH)
                    if len(spp)>0:
                        linestyle=next(linestyles)
                        symbol=next(symbols)
                        if scale=='RHE':
                            plot_stuff(spp,DATA,fit,0,skip,voltage_mode='first',take_log=True,linestyle=linestyle,symbol=symbol,fit_tafel=fit_tafel) #,convert='RHE_to_SHE')
                        elif scale=='SHE':
                            plot_stuff(spp,DATA,fit,0,skip,voltage_mode='first',take_log=True,linestyle=linestyle,symbol=symbol,convert='RHE_TO_SHE',fit_tafel=fit_tafel) #,convert='RHE_to_SHE')


    plt.xlim([-2,1])
    plt.ylim([-6,2])
#    data=np.loadtxt('H2_OD-Cu1_01M_KOH.csv',delimiter=',')
#    plt.plot(data[:,0],np.log10(data[:,1]),'-o',label='H2_OD-Cu1_01M_KOH')
#    data=np.loadtxt('H2_OD-Cu2_01M_KOH.csv',delimiter=',')
#    plt.plot(data[:,0],np.log10(data[:,1]),'-o',label='H2_OD-Cu2_01M_KOH')
#    data=np.loadtxt('H2_np-Cu_01M_KOH.csv',delimiter=',')
#    plt.plot(data[:,0],np.log10(data[:,1]),'-o',label='H2_np-Cu_01M_KOH')
    plt.show()

#plot_data(reference=['hori','jaramillo'],species=['C1-sum','C2-sum','H$_2$'],pH=['6.8','13'],scale='RHE',system=['all'],fit_tafel=True)
#plot_data(reference=['hori','jaramillo','kanan'],species=['C1','HCOO','C2+-sum','HCOO','H$_2$','CO'],pH=['6.8','7.2'],scale='RHE',system=['all'])
plot_data(reference=['hori','jaramillo'],species=['H$_2$'],pH=['6.8','13'],scale='SHE')

sys.exit()

################
#pH effect on HER
################

fig = plt.figure() #1, figsize=figsize)
ax = fig.add_subplot(1,1,1)
voltage_all=[]
current_all=[]

plot_stuff([-1],DATA_jr2,2,0,0,linestyle='-',symbol='o',voltage_mode='first',take_log=True)
plot_stuff([1],DATA_lr,2,0,0,linestyle=':',symbol="D",voltage_mode='first',take_log=True)
plt.xlabel('$U_{\mathrm{RHE}}$ [V] ')
plt.xlim([-2,0])
plt.ylim([-6,1])

plt.show()
sys.exit()


################
#all C1's and CO and H2
################

fig = plt.figure() #1, figsize=figsize)
ax = fig.add_subplot(1,1,1)
voltage_all=[]
current_all=[]
#SHE
#plot_stuff([2,3,4,5],DATA_h2s,2,0,voltage_mode='first',take_log=True,linestyle='-') #,skip=4) #
#plot_stuff([1],DATA_js,2,0,0,linestyle=':',symbol='D')
#plot_stuff([1,2,5],DATA_kr,2,0,0,voltage_mode='first',take_log=True,linestyle='--',symbol='d',convert='RHE_to_SHE')
#plt.xlabel('$U_{\mathrm{SHE}}$ [V] ')
#RHE
plot_stuff([2,3,4,5],DATA_h2s,2,0,voltage_mode='first',take_log=True,linestyle='-',convert='SHE_TO_RHE') #,skip=4) #
#plot_stuff([1],DATA_js,2,0,0,linestyle=':',symbol='D',convert='SHE_to_RHE')
plot_stuff([1,2,5],DATA_kr,2,0,0,voltage_mode='first',take_log=True,linestyle='--',symbol='d') #,convert='RHE_to_SHE')
plot_stuff([1,2,3,4],DATA_jr2,2,0,0,voltage_mode='first',linestyle=':',symbol="D",take_log=True)
#plot_stuff([1,3,4],DATA_ss,2,0,0,voltage_mode='first',linestyle='-',symbol="p",take_log=True,convert='SHE_TO_RHE')
plot_stuff([1,3,4],DATA_ss_01,2,0,0,voltage_mode='first',linestyle='--',symbol="p",take_log=True,convert='SHE_TO_RHE')
plt.xlabel('$U_{\mathrm{RHE}}$ [V] ')
plt.xlim([-2,0])
plt.ylim([-6,1])
#plot_stuff([3,17],DATA_hs,0,0)
#plot_tafel(voltage_all,current_all,'k')
#fig.savefig('jaramillo_current_densities_she_rhe.pdf')

###############
#BUFFER CONC dependence
###############


fig = plt.figure()
ax = fig.add_subplot(1,1,1)
voltage_all=[]
current_all=[]

plot_stuff([1,3,4],DATA_ss_005,0,1,0,voltage_mode='first',linestyle='-',symbol="p",take_log=True,convert='SHE_TO_RHE')
plot_stuff([1,3,4],DATA_ss_01,0,1,0,voltage_mode='first',linestyle='--',symbol="o",take_log=True,convert='SHE_TO_RHE')
plot_stuff([1,3,4],DATA_ss_02,0,1,0,voltage_mode='first',linestyle=':',symbol="d",take_log=True,convert='SHE_TO_RHE')

plt.xlabel('$U_{\mathrm{RHE}}$ [V] ')
plt.xlim([-2,0])
plt.ylim([-6,1])


################
#all C2's
################

fig = plt.figure() #1, figsize=figsize)
ax = fig.add_subplot(1,1,1)
voltage_all=[]
current_all=[]

plot_stuff([3,5],DATA_jr,2,0,0,linestyle=':',symbol="D")
plot_stuff([2],DATA_ss_01,2,0,0,voltage_mode='first',linestyle='--',symbol="p",take_log=True,convert='SHE_TO_RHE')
plot_stuff([1],DATA_h2s,2,0,voltage_mode='first',take_log=True,linestyle='-',convert='SHE_TO_RHE') #,skip=4) #
plot_stuff([3],DATA_kr,2,0,0,voltage_mode='first',take_log=True,linestyle='--',symbol='d') #,convert='RHE_to_SHE')

plt.xlabel('$U_{\mathrm{RHE}}$ [V] ')
plt.xlim([-2,0])
plt.ylim([-6,1])


#all HORI pH = 7
fig = plt.figure() #1, figsize=figsize)
ax = fig.add_subplot(1,1,1)
voltage_all=[]
current_all=[]
plot_stuff(range(1,6),DATA_h2s,2,0,voltage_mode='first',take_log=True) #,skip=4) #
plot_stuff([1],DATA_js,1,1,1)
#plot_stuff([3,17],DATA_hs,0,0)
#plot_tafel(voltage_all,current_all,'k')
#fig.savefig('jaramillo_current_densities_she_rhe.pdf')
plt.show()


#C2's
fig = plt.figure(1, figsize=figsize)
ax = fig.add_subplot(2,2,1)
voltage_all=[]
current_all=[]
plot_stuff([1,3,5,7,9,11,13],DATA_hs,0,0)
plot_tafel(voltage_all,current_all,'k')
plt.xlabel('$U_{\mathrm{SHE}}$ [V] ')
plt.ylim((-3,1))
plt.xlim((-1.4,-1.05))


ax2 = fig.add_subplot(2,2,2)
# voltage_all=[]
# current_all=[]
plot_stuff([1,3,5,7,9,11,13],DATA_hr,1,0)
#plot_tafel(voltage_all,current_all,'k')
plt.xlabel('$U_{\mathrm{RHE}}$ [V] ')
plt.ylim((-3,1))
plt.xlim((-0.9,-0.6))


#C1's
ax3=fig.add_subplot(2,2,3)
plot_stuff([15,17,19,21,23,25,27],DATA_hs,1,0)
plt.xlabel('$U_{\mathrm{SHE}}$ [V] ')
plt.ylim((-3,1))
plt.xlim((-1.4,-1.05))



ax4=fig.add_subplot(2,2,4)
plot_stuff([15,17,19,21,23,25,27],DATA_hr,1,0)
plt.xlabel('$U_{\mathrm{RHE}}$ [V] ')
plt.ylim((-3,1))
plt.xlim((-0.9,-0.6))

fig.savefig('hori_current_densities_she_rhe.pdf')


#Jaramillo  
#C2's
fig1 = plt.figure(2, figsize=figsize)
ax=fig1.add_subplot(2,2,1)
plot_stuff([3,5,9,11,13],DATA_js,1,1,1)
plt.xlabel('$U_{\mathrm{SHE}}$ [V] ')
plt.ylim((-4.5,2))
plt.xlim((-1.7,-0.9))


ax2=fig1.add_subplot(2,2,2)
plot_stuff([3,5,9,11,13],DATA_jr,1,1,1)

plt.xlabel('$U_{\mathrm{RHE}}$ [V] ')
plt.ylim((-4.5,2))
plt.xlim((-1.2,-0.2))


ax=fig1.add_subplot(2,2,3)
plot_stuff([1,7],DATA_js,1,1,1)

plt.xlabel('$U_{\mathrm{SHE}}$ [V] ')
plt.ylim((-4.5,2))
plt.xlim((-1.7,-0.9))


#C1's
ax2=fig1.add_subplot(2,2,4)
plot_stuff([1,7],DATA_jr,1,1,1)
plt.xlabel('$U_{\mathrm{RHE}}$ [V] ')
plt.ylim((-4.5,2))
plt.xlim((-1.2,-0.2))


fig1.savefig('jaramillo_current_densities_she_rhe.pdf')



#Jaramillo & Hori   
#C2's
fig2 = plt.figure(3, figsize=figsize)
ax=fig2.add_subplot(2,2,1)
plot_stuff([3,5,9,11,13],DATA_js,1,1,1)
voltage_all=[]
current_all=[]
plot_stuff([1,3,5,7,9,11,13],DATA_hs,0,0)
plot_tafel(voltage_all,current_all,'k')
plt.xlabel('$U_{\mathrm{SHE}}$ [V] ')
plt.ylim((-4.5,2))
plt.xlim((-1.7,-0.8))


ax2=fig2.add_subplot(2,2,2)
plot_stuff([3,5,9,11,13],DATA_jr,1,1,1)
plot_stuff([1,3,5,7,9,11,13],DATA_hr,1,0)
plt.xlabel('$U_{\mathrm{RHE}}$ [V] ')
plt.ylim((-4.5,2))
plt.xlim((-1.2,-0.2))


ax=fig2.add_subplot(2,2,3)
plot_stuff([1,7],DATA_js,1,1,1)
plot_stuff([15,17,19,21,23,25,27],DATA_hs,1,0)

plt.xlabel('$U_{\mathrm{SHE}}$ [V] ')
plt.ylim((-4.5,2))
plt.xlim((-1.7,-0.8))


#C1's
ax2=fig2.add_subplot(2,2,4)
plot_stuff([1,7],DATA_jr,1,1,1)
plot_stuff([15,17,19,21,23,25,27],DATA_hr,1,0)
plt.xlabel('$U_{\mathrm{RHE}}$ [V] ')
plt.ylim((-4.5,2))
plt.xlim((-1.2,-0.2))
    
fig2.savefig('jaramillo_hori_current_densities_she_rhe.pdf')







