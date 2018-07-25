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
import os


class Object(object):
    pass

class EXPDATA():
    def __init__(self):

        font = {'family' : 'serif',
                'serif':['Times'],
                 'size'   : 16}
        rc('font', **font)
        rc('text', usetex=True)   
        rc('xtick', labelsize=16)
        rc('ytick', labelsize=16)

        X=7
        figsize = (3*X,2*X)
        
        self.color = ['r','darkorange','gold','g','c','b','violet','r','darkorange','gold','g','c','b','violet']
        linestyle = [':',':','-','-','-','-','-','-','-','-','-','-','-','-','-']
        self.marker=['o','o','o','o','o','o','o','v','v','v','v','v','v','v']
        
        root=os.getcwd()
        os.chdir('/'.join(os.path.realpath(__file__).split('/')[:-1])+'/../data/CO2R')

        self.DATA_hs=self.get_data('COR_hori_ph_normalized_SHE.csv')
        self.DATA_hr=self.get_data('COR_hori_ph_normalized_RHE.csv') #hori_data_RHE.csv')

        self.DATA_jr=self.get_data('CO2R_jaramillo_normalized_rhe.csv')

        self.DATA_jr_NF=self.get_data('CO2R_jaramillo_nfCu_normalized_RHE.csv') #nanoflowers

        self.DATA_jr_NC=self.get_data('CO2R_CuCubes.csv') #nanocubes

        self.DATA_lr_HER=self.get_data('COR_lei_HER_data.csv')
        self.DATA_lr_NF=self.get_data('COR_lei_NF_Cu.csv')


        self.DATA_h2s=self.get_data('CO2R_hori_normalized_SHE.csv')
       
        self.DATA_wu=self.get_data('CO2R_wuttig_normalized_SHE.csv')
        self.DATA_du=self.get_data('CO2R_dunwell_normalized_RHE.csv')

        self.DATA_kr=self.get_data('CO2R_kanan_normalized_RHE.csv')
        self.DATA_wr_NW=self.get_data('COR_Wang_NW.csv')
        self.DATA_wr_NW_all=self.get_data('COR_Wang_NW_all.csv')
        
        self.DATA_ss_005=self.get_data('strasser_jpart_from_prate_KHCO3_005_CO2R_SHE.csv')
        self.DATA_ss_02=self.get_data('strasser_jpart_from_prate_KHCO3_02_CO2R_SHE.csv')
        self.DATA_ss_01=self.get_data('strasser_jpart_from_prate_KHCO3_01_CO2R_SHE.csv')
        self.DATA_ss_005_FE=self.get_data('strasser_jpart_KHCO3_005_CO2R_NHE.csv')
        self.DATA_ss_01_FE=self.get_data('strasser_jpart_KHCO3_01_CO2R_NHE.csv')
        self.DATA_ss_02_FE=self.get_data('strasser_jpart_KHCO3_02_CO2R_NHE.csv')
        
        self.DATA_lr=self.get_data('COR_lei_normalized_rhe.csv')

        self.DATA_cas=self.get_data('CO2R_Au_carlos_SHE.csv')
        self.DATA_cas2=self.get_data('CO2R_Au_carlos_SHE_pH3.csv')
#        self.DATA_lr2=self.get_data('COR_lei_high_surface_x380_rhe.csv')
        os.chdir(root)
        


    def get_data(self,filename):
        o = open('CSV/'+filename,'rU')
        DATA = Object()
        DATA.data = np.array(list(csv.reader(o)))
        print filename
        print 'check',np.shape(DATA.data)
        DATA.label=[str(d) for d in DATA.data[0,:]] #data[0,1:4]
        DATA.data = np.delete(DATA.data,0,0) # delete header
        return DATA

    def plot_stuff(self,list_of_data,DATA,fit,joinlines, skip={}, ax=None,voltage_mode='previous',take_log=False,\
            linestyle=None,symbol=None,convert=None,fit_tafel=False,legend=False,msize=12,color=None):
        #if voltage_mode='previous' take previous point before list_of_data as voltages, if 'first', take first column for all data
        #if fit=0 Tafel, fit=1 poly, fit=2 poly no points
        #fit_tafel: New tafel fitting based on error estimation due to curvature
        n=0
        if ax is not None:
            bf=ax
        else:
            bf=plt
#        if take_log:
#            func=bf.semilogy
#        else:
#            func=bf.plot
        func=bf.semilogy

        symbols=self.symbols
        for j in list_of_data: #[1,3,5,7,9,11,13]:
            if symbol is not None: #and j==0:
                cmarker=symbol
            #elif symbol is not None:
            #    cmarker=next(symbols)
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
#            if any(current)<0:
#                current=[-c for c in current]
            if not take_log:
                current=[10**c for c in current]
            for k in range(0,len(voltage)):
                voltage_all.append(voltage[k])
                #if take_log:
                #    current[k]=np.log10(current[k])
                current_all.append(current[k])
            linestyle_2 = "None"
            if joinlines:
                linestyle_2=":"
            if linestyle is None:
                linestyle = "None"
                if joinlines:
                    linestyle=":"
            if color is None:
                self.color[n]=self.get_color(DATA.label[j])
            else:
                self.color[n]=color
            val=re.findall('pH[ ]*=[ ]*(\d+)',DATA.label[j])
            if len(val)>0:
                pH=float(val[0])
            else:
                pH=7.0
            if symbol is None:
                cmarker=self.marker[n]
            if convert is not None:
                if convert.split('_')[0]=='RHE':
                    shift=-0.0592*pH
                else:
                    shift=+0.0592*pH
            else:
                shift=0.0
            voltage=[v+shift for v in voltage]
    
            #label
            #labels.append(DATA.label[j])
            #legend_labels = numpy.concatenate([label_row_1, label_j_1, label_empty * 3, label_j_2, label_empty * 3, label_j_3, label_empty * 3])
    
            #if fit_tafel:
            #    tafel_fitting(voltage,current)
            if not fit>2:
#                plt.figure()
                func(voltage, current, color=self.color[n], linestyle=linestyle_2, marker = cmarker, markersize=msize,label=DATA.label[j])  #linestyle = ':',
#                plt.show()
#                sys.exit()
            if fit_tafel:
                fit=1
                linestyle=':'
            if fit>0:
                if len(voltage)>1:
                    species=DATA.label[j].split(',')[0].strip()
                    skip_val=0
                    for sk in skip:
                        if sk in species:
                            skip_val=skip[sk]
                    if fit==1:
                        self.plot_tafel(bf,voltage,current,self.color[n],skip_val,linestyle=linestyle,take_log=take_log)
                    if fit>1:
                        V=voltage
                        J=current
                        linecol=self.color[n]
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
                            func(xfine,p(xfine),linestyle=linestyle,color=linecol,label=DATA.label[j])
                        else:
                            print 'plotting fit',DATA.label[j],linestyle
                            print 'data',zip(V,J,p(V))
                            func(xfine,p(xfine),linestyle=linestyle,color=linecol)
                            
            n=n+1   
        #plt.ylabel('log$j$ [mA/cm$^2_{\mathrm{real}}$] ')
        #plt.ylim((-3,1))
        #plt.xlim((-1.4,-1.05))
        if legend:
            leg = ax.legend(loc=1, 
                  ncol=1, fontsize=8, numpoints=1) #, fancybox=True, shadow=False, loc='upper center', bbox_to_anchor=(0.5, 1.05))
            for line,text in zip(leg.get_lines(), leg.get_texts()):
                text.set_color(line.get_color())
    
    def plot_tafel(self,ax,V,J,linecol,skip=0,linestyle='-',take_log=True):
        # fit lines
        if take_log:
            J=np.log10(J)
        data=[V,J]
        data=map(list,zip(*data))
        data=np.array(data)
        #print 'before',data
        #data=sorted(data,key=lambda l:l[1], reverse=True)
        #data=np.array(data)
        #print 'after',data
        data.view(data.dtype.str+','+data.dtype.str).sort(order=['f0'], axis=0)
        V=data[:,0]
        J=data[:,1]
        if skip:
            V=V[skip:]
            J=J[skip:]
        #A = array([V, ones(len(V))])
        #w = linalg.lstsq(A.T,J )[0]
        #line_tafel = w[0]*np.array(V)+w[1]
        #tafel_slope = 1000/w[0]
        p=np.polyfit(V,J,1)
        z=np.poly1d(p)
        tafel_slope=-1000./p[1]
        if 300<tafel_slope<500:
            return ax
        ax.semilogy(V,10**z(V),color=linecol, linestyle=linestyle)
        ax.text(V[0],10**z(V[0]),str(int(tafel_slope))+'mV/dec', color=linecol, fontsize=11)  #
        return ax
    
    def plot_data(self,ax=None,reference=['all'],species=['all'],pH=['all'],ci_bic=['all'],scale='RHE',reaction='all',\
            system=['all'],coloring='species',fit_tafel=False,only_points=False,\
            take_log=True,marker=None,legend=False,msize=None,color=None):
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
        only_points: if only points without tafel or fit should be plotted
        take_log: if true this plots the logarithm of the data
        legend: if true, plot legend
        msize: force size of markers
        color: force color of plot (otherwise colored by product)
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
        species=['C$_2$' if sp=='C2-sum' else sp for sp in species]
        species=['C$_{2+}$' if sp=='C2+-sum' else sp for sp in species]
    
        print('Plotting species = {}'.format(species))
    
#        plt.figure()
    
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
                        if h.startswith(sp) and anyph(pH,h) and reac(h) and system(h) and i not in heads:
                            heads+=[i]
            return heads
    
        linestyles=cycle(['-','-',':','-.'])
        self.symbols=cycle(['d','D','o','x','p','*','v','h','1','2'])
        symbols=self.symbols
    
        if fit_tafel:
            fit=1
        else:
            fit=2

        if only_points:
            fit=0
        
    
    #    data_labels=[','.join(a.label[1].split(',')[1:]) for a in [self.DATA_ss_01,self.DATA_ss_005,self.DATA_ss_02,self.DATA_ss_005_FE,self.DATA_ss_01_FE,self.DATA_ss_02_FE,self.DATA_jr,self.DATA_kr,self.DATA_lr,self.DATA_h2s,self.DATA_hr,self.DATA_hs,self.DATA_]]
        data_labels=[','.join(a.label[1].split(',')[1:]) for a in [self.DATA_ss_01,self.DATA_ss_005,self.DATA_ss_02,self.DATA_ss_005_FE,self.DATA_ss_01_FE,self.DATA_ss_02_FE,self.DATA_jr,self.DATA_jr_NF,self.DATA_jr_NC,self.DATA_kr,self.DATA_lr,self.DATA_lr_NF,self.DATA_lr_HER,self.DATA_h2s,self.DATA_hr,self.DATA_hs,self.DATA_wr_NW,self.DATA_wr_NW_all,self.DATA_wu,self.DATA_du,self.DATA_cas,self.DATA_cas2]]
    
    
    
    
        for cpH in pH:
            skip_dict={}
            for label in data_labels:
                if 'Hori' in label:
                    #skip_dict[label]={'HCOO':   5,\
                    #                  'CO':     9,\
                    #                  'H$_2$':  10}
                    skip_dict[label]={'HCOO':   0,\
                                      'CO':     0,\
                                      'H$_2$':  0}
                elif 'Wang' in label:
                    skip_dict[label]={'C$_2$H$_4$': 0}
                elif 'Jaramillo' in label:
                    #skip_dict[label]={'HCOO':   6,\
                    #                  'CH$_4$': 2,\
                    #                  'H$_2$':  5,\
                    #                  'C$_2$H$_4$':3,\
                    #                  'EtOH':3,\
                    #                  'Acetate':0,\
                    #                  'propol':0}
                    skip_dict[label]={'HCOO':   0,\
                                      'CH$_4$': 0,\
                                      'H$_2$':  0,\
                                      'C$_2$H$_4$':0,\
                                      'EtOH':0,\
                                      'Acetate':0,\
                                      'propol':0}
                    if 'pc-Au' in label and 'pH = 3.0' in label:
                        skip_dict[label]={'CO': 3}
                    elif 'pc-Au' in label and 'pH = 6.8' in label:
                        skip_dict[label]={'CO': 0}
                elif 'Kanan' in label:
                    skip_dict[label]={'HCOO':   9,\
                                      'CO':     10,\
                                      'H$_2$':     9,\
                                      'C$_2$H$_4$':2,\
                                      'C$_2$H$_6$':3}
                elif 'Wuttig' in label or 'Dunwell' in label:
                    skip_dict[label]={}
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
                        DATA=self.DATA_ss_01
                        name=','.join(DATA.label[1].split(',')[1:])
                        skip=skip_dict[name]
                        spp=s2i(species,DATA,pH=cpH)
                        if len(spp)>0:
                            linestyle=next(linestyles)
                            if marker is None:
                                symbol=next(symbols)
                            else:
                                symbol=marker
                            if scale=='RHE':
                                self.plot_stuff(spp,DATA,fit,0,skip,ax=ax,linestyle=linestyle,symbol=symbol,\
                                    take_log=take_log,convert='SHE_TO_RHE',fit_tafel=fit_tafel,legend=legend,\
                                    msize=msize,color=color)
                            elif scale=='SHE':
                                self.plot_stuff(spp,DATA,fit,0,skip,ax=ax,linestyle=linestyle,symbol=symbol,\
                                    take_log=take_log,fit_tafel=fit_tafel,legend=legend,\
                                    msize=msize,color=color)
                        ##alternatively from FE's
                        #DATA=self.DATA_ss_01_FE
                        #name=','.join(DATA.label[1].split(',')[1:])
                        #skip=skip_dict[name]
                        #spp=s2i(species,DATA,pH=cpH)
                        #if len(spp)>0:
                        #    linestyle=next(linestyles)
                        #    if marker is None:
                        #        symbol=next(symbols)
                        #    else:
                        #        symbol=marker
                        #    if scale=='RHE':
                        #        self.plot_stuff(spp,DATA,fit,0,skip,ax=ax,linestyle=linestyle,symbol=symbol,\
                        #            voltage_mode='first',take_log=take_log,convert='SHE_TO_RHE',fit_tafel=fit_tafel,\
                        #            msize=msize,legend=legend,color=color)
                        #    else:
                        #        self.plot_stuff(spp,DATA,fit,0,skip,ax=ax,linestyle=linestyle,symbol=symbol,\
                        #            voltage_mode='first',take_log=take_log,fit_tafel=fit_tafel,legend=legend,\
                        #            msize=msize,color=color)
                    if bic in ['0.05','all']:
                        DATA=self.DATA_ss_005
                        name=','.join(DATA.label[1].split(',')[1:])
                        skip=skip_dict[name]
                        spp=s2i(species,DATA,pH=cpH)
                        if len(spp)>0:
                            linestyle=next(linestyles)
                            if marker is None:
                                symbol=next(symbols)
                            else:
                                symbol=marker
                            if scale=='RHE':
                                self.plot_stuff(spp,DATA,fit,0,skip,ax=ax,linestyle=linestyle,symbol=symbol,\
                                    take_log=take_log,convert='SHE_TO_RHE',fit_tafel=fit_tafel,legend=legend,\
                                    msize=msize,color=color)
                            elif scale=='SHE':
                                self.plot_stuff(spp,DATA,fit,0,skip,ax=ax,linestyle=linestyle,symbol=symbol,\
                                    take_log=take_log,fit_tafel=fit_tafel,legend=legend,msize=msize,color=color)
                        ##alternatively from FE's
                        #DATA=self.DATA_ss_005_FE
                        #name=','.join(DATA.label[1].split(',')[1:])
                        #skip=skip_dict[name]
                        #spp=s2i(species,DATA,pH=cpH)
                        #if len(spp)>0:
                        #    linestyle=next(linestyles)
                        #    if marker is None:
                        #        symbol=next(symbols)
                        #    else:
                        #        symbol=marker
                        #    if scale=='RHE':
                        #        self.plot_stuff(spp,DATA,fit,0,skip,ax=ax,linestyle=linestyle,symbol=symbol,\
                        #            take_log=take_log,voltage_mode='first',convert='SHE_TO_RHE',fit_tafel=fit_tafel,\
                        #            legend=legend,msize=msize,color=color)
                        #    else:
                        #        self.plot_stuff(spp,DATA,fit,0,skip,ax=ax,linestyle=linestyle,symbol=symbol,\
                        #            voltage_mode='first',take_log=take_log,fit_tafel=fit_tafel,legend=legend,\
                        #            msize=msize,color=color)
                    if bic in ['0.2','all']:
                        DATA=self.DATA_ss_02
                        name=','.join(DATA.label[1].split(',')[1:])
                        skip=skip_dict[name]
                        spp=s2i(species,DATA,pH=cpH)
                        if len(spp)>0:
                           # sys.exit()
                            linestyle=next(linestyles)
                            if marker is None:
                                symbol=next(symbols)
                            else:
                                symbol=marker
                            if scale=='RHE':
                                self.plot_stuff(spp,DATA,fit,0,skip,ax=ax,linestyle=linestyle,symbol=symbol,\
                                    take_log=take_log,convert='SHE_TO_RHE',fit_tafel=fit_tafel,legend=legend,\
                                    msize=msize,color=color)
                            elif scale=='SHE':
                                self.plot_stuff(spp,DATA,fit,0,skip,ax=ax,linestyle=linestyle,symbol=symbol,\
                                    take_log=take_log,fit_tafel=fit_tafel,legend=legend,msize=msize,color=color)
                        ##alternatively from FE's
                        #DATA=self.DATA_ss_02_FE
                        #name=','.join(DATA.label[1].split(',')[1:])
                        #skip=skip_dict[name]
                        #spp=s2i(species,DATA,pH=cpH)
                        #if len(spp)>0:
                        #    linestyle=next(linestyles)
                        #    if marker is None:
                        #        symbol=next(symbols)
                        #    else:
                        #        symbol=marker
                        #    if scale=='RHE':
                        #        self.plot_stuff(spp,DATA,fit,0,skip,ax=ax,linestyle=linestyle,symbol=symbol,\
                        #            take_log=take_log,voltage_mode='first',convert='SHE_TO_RHE',fit_tafel=fit_tafel,\
                        #            legend=legend,msize=msize,color=color)
                        #    else:
                        #        self.plot_stuff(spp,DATA,fit,0,skip,ax=ax,linestyle=linestyle,symbol=symbol,\
                        #            take_log=take_log,voltage_mode='first',fit_tafel=fit_tafel,legend=legend,\
                        #            msize=msize,color=color)
                if bic in ['0.1','all']:
                    if isref('jaramillo-her'):
                        DATA=self.DATA_lr_HER
                        name=','.join(DATA.label[1].split(',')[1:])
                        skip=skip_dict[name]
                        spp=s2i(species,DATA,pH=cpH)
                        if len(spp)>0:
                            linestyle=next(linestyles)
                            if marker is None:
                                symbol=next(symbols)
                            else:
                                symbol=marker
                            if scale=='RHE':
                                self.plot_stuff(spp,DATA,fit,0,skip,ax=ax,linestyle=linestyle,symbol=symbol,\
                                    voltage_mode='previous',take_log=take_log,fit_tafel=fit_tafel,legend=legend,\
                                    msize=msize,color=color)
                            elif scale=='SHE':
                                self.plot_stuff(spp,DATA,fit,0,skip,ax=ax,linestyle=linestyle,symbol=symbol,\
                                    voltage_mode='previous',take_log=take_log,convert='RHE_TO_SHE',fit_tafel=fit_tafel,\
                                    msize=msize,legend=legend,color=color)
                    if isref('jaramillo'):
                        DATA=self.DATA_jr
                        name=','.join(DATA.label[1].split(',')[1:])
                        skip=skip_dict[name]
                        spp=s2i(species,DATA,pH=cpH)
                        if len(spp)>0:
                            linestyle=next(linestyles)
                            if marker is None:
                                symbol=next(symbols)
                            else:
                                symbol=marker
                            if scale=='RHE':
                                self.plot_stuff(spp,DATA,fit,0,skip,ax=ax,linestyle=linestyle,symbol=symbol,\
                                    voltage_mode='first',take_log=take_log,fit_tafel=fit_tafel,legend=legend,\
                                    msize=msize,color=color)
                            elif scale=='SHE':
                                self.plot_stuff(spp,DATA,fit,0,skip,ax=ax,linestyle=linestyle,symbol=symbol,\
                                    voltage_mode='first',take_log=take_log,convert='RHE_TO_SHE',fit_tafel=fit_tafel,\
                                    msize=msize,legend=legend,color=color)
                        DATA=self.DATA_lr
                        name=','.join(DATA.label[1].split(',')[1:])
                        skip=skip_dict[name]
                        spp=s2i(species,DATA,pH=cpH)
                        if len(spp)>0:
                            linestyle=next(linestyles)
                            if marker is None:
                                symbol=next(symbols)
                            else:
                                symbol=marker
                            if scale=='RHE':
                                self.plot_stuff(spp,DATA,fit,0,skip,ax=ax,linestyle=linestyle,symbol=symbol,\
                                    voltage_mode='first',take_log=take_log,fit_tafel=fit_tafel,legend=legend,\
                                    msize=msize,color=color)
                            else:
                                self.plot_stuff(spp,DATA,fit,0,skip,ax=ax,linestyle=linestyle,symbol=symbol,\
                                    voltage_mode='first',take_log=take_log,convert='RHE_TO_SHE',fit_tafel=fit_tafel,\
                                    msize=msize,legend=legend,color=color)
                        DATA=self.DATA_lr_NF                                                                         ##
                        name=','.join(DATA.label[1].split(',')[1:])
                        skip=skip_dict[name]
                        spp=s2i(species,DATA,pH=cpH)
                        if len(spp)>0:
                            linestyle=next(linestyles)
                            symbol=next(symbols)
                            if scale=='RHE':
                                self.plot_stuff(spp,DATA,fit,0,skip,linestyle=linestyle,symbol=symbol,\
                                    voltage_mode='first',take_log=True,fit_tafel=fit_tafel,legend=legend,ax=ax,msize=msize,color=color)
                            else:
                                self.plot_stuff(spp,DATA,fit,0,skip,linestyle=linestyle,symbol=symbol,\
                                    voltage_mode='first',take_log=True,convert='RHE_TO_SHE',fit_tafel=fit_tafel,legend=legend,\
                                    msize=msize,ax=ax,color=color)
                        DATA=self.DATA_jr_NF
                        name=','.join(DATA.label[1].split(',')[1:])
                        skip=skip_dict[name]
                        spp=s2i(species,DATA,pH=cpH)
                        if len(spp)>0:
                            linestyle=next(linestyles)
                            symbol=next(symbols)
                            if scale=='RHE':
                                self.plot_stuff(spp,DATA,fit,0,skip,linestyle=linestyle,symbol=symbol,\
                                    voltage_mode='first',take_log=True,fit_tafel=fit_tafel,legend=legend,ax=ax,msize=msize,color=color)
                            elif scale=='SHE':
                                self.plot_stuff(spp,DATA,fit,0,skip,linestyle=linestyle,symbol=symbol,\
                                    voltage_mode='first',take_log=True,convert='RHE_TO_SHE',fit_tafel=fit_tafel,legend=legend,\
                                    msize=msize,ax=ax,color=color)
                        DATA=self.DATA_jr_NC
                        name=','.join(DATA.label[1].split(',')[1:])
                        skip=skip_dict[name]
                        spp=s2i(species,DATA,pH=cpH)
                        if len(spp)>0:
                            linestyle=next(linestyles)
                            symbol=next(symbols)
                            if scale=='RHE':
                                self.plot_stuff(spp,DATA,fit,0,skip,linestyle=linestyle,symbol=symbol,\
                                    voltage_mode='first',take_log=True,fit_tafel=fit_tafel,legend=legend,ax=ax,msize=msize,color=color)
                            elif scale=='SHE':
                                self.plot_stuff(spp,DATA,fit,0,skip,linestyle=linestyle,symbol=symbol,\
                                    voltage_mode='first',take_log=True,convert='RHE_TO_SHE',fit_tafel=fit_tafel,legend=legend,\
                                    msize=msize,ax=ax,color=color)
                        DATA=self.DATA_cas
                        name=','.join(DATA.label[1].split(',')[1:])
                        skip=skip_dict[name]
                        spp=s2i(species,DATA,pH=cpH)
                        if len(spp)>0:
                            linestyle=next(linestyles)
                            symbol=next(symbols)
                            if scale=='SHE':
                                self.plot_stuff(spp,DATA,fit,0,skip,linestyle=linestyle,symbol=symbol,\
                                    voltage_mode='first',take_log=True,fit_tafel=fit_tafel,legend=legend,ax=ax,msize=msize,color=color)
                            elif scale=='RHE':
                                self.plot_stuff(spp,DATA,fit,0,skip,linestyle=linestyle,symbol=symbol,\
                                    voltage_mode='first',take_log=True,convert='SHE_TO_RHE',fit_tafel=fit_tafel,legend=legend,\
                                    msize=msize,ax=ax,color=color)
                        DATA=self.DATA_cas2
                        name=','.join(DATA.label[1].split(',')[1:])
                        skip=skip_dict[name]
                        spp=s2i(species,DATA,pH=cpH)
                        if len(spp)>0:
                            linestyle=next(linestyles)
                            symbol=next(symbols)
                            if scale=='SHE':
                                self.plot_stuff(spp,DATA,fit,0,skip,linestyle=linestyle,symbol=symbol,\
                                    voltage_mode='first',take_log=True,fit_tafel=fit_tafel,legend=legend,ax=ax,msize=msize,color=color)
                            elif scale=='RHE':
                                self.plot_stuff(spp,DATA,fit,0,skip,linestyle=linestyle,symbol=symbol,\
                                    voltage_mode='first',take_log=True,convert='SHE_TO_RHE',fit_tafel=fit_tafel,legend=legend,\
                                    msize=msize,ax=ax,color=color)
#                        DATA=self.DATA_lr2
#                        name=','.join(DATA.label[1].split(',')[1:])
#                        skip=skip_dict[name]
#                        spp=s2i(species,DATA,pH=cpH)
#                        if len(spp)>0:
#                            linestyle=next(linestyles)
#                            symbol=next(symbols)
#                            if scale=='RHE':
#                                self.plot_stuff(spp,DATA,fit,0,skip,ax=ax,linestyle=linestyle,symbol=symbol,voltage_mode='first',take_log=take_log,fit_tafel=fit_tafel)
#                            else:
#                                self.plot_stuff(spp,DATA,fit,0,skip,ax=ax,linestyle=linestyle,symbol=symbol,voltage_mode='first',take_log=take_log,convert='RHE_TO_SHE',fit_tafel=fit_tafel)
                    if isref('wang'):
                        DATA=self.DATA_wr_NW_all #DATA_wr_NW
                        name=','.join(DATA.label[1].split(',')[1:])
                        skip=skip_dict[name]
                        spp=s2i(species,DATA,pH=cpH)
                        if len(spp)>0:
                            linestyle=next(linestyles)
                            symbol=next(symbols)
                            if scale=='RHE':
                                plot_stuff(spp,DATA,fit,0,skip,linestyle=linestyle,symbol=symbol,\
                                    voltage_mode='previous',take_log=True,fit_tafel=fit_tafel,legend=legend,ax=ax,msize=msize,color=color)
                            elif scale=='SHE':
                                plot_stuff(spp,DATA,fit,0,skip,linestyle=linestyle,symbol=symbol,\
                                    voltage_mode='previous',take_log=True,convert='RHE_TO_SHE',fit_tafel=fit_tafel,legend=legend,ax=ax,msize=msize,color=color)
                    if isref('hori'):
                        DATA=self.DATA_h2s
                        name=','.join(DATA.label[1].split(',')[1:])
                        skip=skip_dict[name]
                        spp=s2i(species,DATA,pH=cpH)
                        if len(spp)>0:
                            linestyle=next(linestyles)
                            if marker is None:
                                symbol=next(symbols)
                            else:
                                symbol=marker
                            if scale=='RHE':
                                self.plot_stuff(spp,DATA,fit,0,skip,ax=ax,voltage_mode='first',take_log=take_log,linestyle=linestyle,symbol=symbol,\
                                    convert='SHE_TO_RHE',fit_tafel=fit_tafel,legend=legend,msize=msize,color=color) #,skip=4) #
                            elif scale=='SHE':
                                self.plot_stuff(spp,DATA,fit,0,skip,ax=ax,voltage_mode='first',take_log=take_log,linestyle=linestyle,symbol=symbol,\
                                    fit_tafel=fit_tafel,legend=legend,msize=msize,color=color) #,skip=4) #
    
                        if scale=='RHE':
                            DATA=self.DATA_hr
                            name=','.join(DATA.label[1].split(',')[1:])
                            skip=skip_dict[name]
                        elif scale=='SHE':
                            DATA=self.DATA_hs
                            name=','.join(DATA.label[1].split(',')[1:])
                            skip=skip_dict[name]
                        spp=s2i(species,DATA,pH=cpH)
                        if len(spp)>0:
                            linestyle=next(linestyles)
                            if marker is None:
                                symbol=next(symbols)
                            else:
                                symbol=marker
                            self.plot_stuff(spp,DATA,fit,0,skip,ax=ax,linestyle=linestyle,symbol=symbol,fit_tafel=fit_tafel,legend=legend,msize=msize,color=color) 
    
                    if isref('wuttig'):
                        DATA=self.DATA_wu
                        name=','.join(DATA.label[1].split(',')[1:])
                        skip=skip_dict[name]
                        spp=s2i(species,DATA,pH=cpH)
                        if len(spp)>0:
                            linestyle=next(linestyles)
                            if marker is None:
                                symbol=next(symbols)
                            else:
                                symbol=marker
                            if scale=='RHE':
                                self.plot_stuff(spp,DATA,fit,0,skip,ax=ax,voltage_mode='first',take_log=False,linestyle=linestyle,\
                                    symbol=symbol,convert='SHE_TO_RHE',fit_tafel=fit_tafel,legend=legend,msize=msize,color=color)
                            elif scale=='SHE':
                                self.plot_stuff(spp,DATA,fit,0,skip,ax=ax,voltage_mode='first',take_log=False,linestyle=linestyle,\
                                    symbol=symbol,fit_tafel=fit_tafel,legend=legend,msize=msize,color=color)

                    if isref('dunwell'):
                        DATA=self.DATA_du
                        name=','.join(DATA.label[1].split(',')[1:])
                        skip=skip_dict[name]
                        spp=s2i(species,DATA,pH=cpH)
                        if len(spp)>0:
                            linestyle=next(linestyles)
                            if marker is None:
                                symbol=next(symbols)
                            else:
                                symbol=marker
                            if scale=='RHE':
                                self.plot_stuff(spp,DATA,fit,0,skip,ax=ax,voltage_mode='first',take_log=False,linestyle=linestyle,\
                                    symbol=symbol,fit_tafel=fit_tafel,legend=legend,msize=msize,color=color)
                            elif scale=='SHE':
                                self.plot_stuff(spp,DATA,fit,0,skip,ax=ax,voltage_mode='first',take_log=False,linestyle=linestyle,\
                                    symbol=symbol,convert='RHE_TO_SHE',fit_tafel=fit_tafel,legend=legend,msize=msize,color=color)
                        

                    if isref('kanan'):
                        DATA=self.DATA_kr
                        name=','.join(DATA.label[1].split(',')[1:])
                        skip=skip_dict[name]
                        spp=s2i(species,DATA,pH=cpH)
                        if len(spp)>0:
                            linestyle=next(linestyles)
                            if marker is None:
                                symbol=next(symbols)
                            else:
                                symbol=marker
                            if scale=='RHE':
                                self.plot_stuff(spp,DATA,fit,0,skip,ax=ax,voltage_mode='first',take_log=take_log,linestyle=linestyle,\
                                    symbol=symbol,fit_tafel=fit_tafel,legend=legend,msize=msize,color=color)
                            elif scale=='SHE':
                                self.plot_stuff(spp,DATA,fit,0,skip,ax=ax,voltage_mode='first',take_log=take_log,linestyle=linestyle,\
                                    symbol=symbol,convert='RHE_TO_SHE',fit_tafel=fit_tafel,legend=legend,msize=msize,color=color)
    
    
#        plt.xlim([-2,1])
#        plt.ylim([-6,2])
#        plt.show()
    def get_color(self,species):
        species=species.strip()
        color=None #'k'
        if any([species.startswith(a) for a in ['n-PrOH','C$_{2+}$']]):
            color='lightblue'
        elif any([species.startswith(a) for a in ['CH3CH2OH','C$_2$H$_4$','CH$_3$COO','etol','C2','EtOH','C$_{2+}$','C$_2$']]):
            color='r'
        elif any([species.startswith(a) for a in ['CH$_4$','C$_1$','CH4','C1']]):
            color='orange'
        elif any([species.startswith(a) for a in ['HCOO','HCOO-','HCOOH']]):
            color='b'
        elif any([species.startswith(a) for a in ['H$_2','H2']]):
            color='k'
        elif any([species.startswith(a) for a in ['CO2','CO$_2$']]):
            color='darkred'
        elif species.startswith('COOH'):
            color='g'
        elif species.startswith('CO'):
            color='olive'
        else:
            pass
        return color
                #color[n]='dark red'
