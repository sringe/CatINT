#!/usr/bin/env python


from numpy import arange,array,ones,linalg
from matplotlib.pyplot import plot,show
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import rc
import re
from itertools import cycle
from .units import *
from scipy.interpolate import UnivariateSpline
import math

import sys
import csv
import sys
import os


class Object(object):
    pass

class EXPDATA():
    def __init__(self):
        itafel=0
        font = {'family' : 'serif',
                'serif':['Times'],
                 'size'   : 16}
        #rc('font', **font)
        #rc('text', usetex=True)   
        #rc('xtick', labelsize=16)
        #rc('ytick', labelsize=16)

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


        #hori on Cu
        self.DATA_h2s=self.get_data('CO2R_hori_normalized_SHE.csv')
        #hori on Au
        self.DATA_hau=self.get_data('CO2R_Au_hori_SHE.csv')
       
        self.DATA_wu=self.get_data('CO2R_wuttig_normalized_SHE.csv')
        self.DATA_du=self.get_data('CO2R_dunwell_normalized_RHE.csv')
        self.DATA_du2=self.get_data('CO2R_dunwell_normalized_RHE_2.csv')

        #kanan on Cu
        self.DATA_kr=self.get_data('CO2R_kanan_normalized_RHE.csv')
        #kanan on Au
        self.DATA_kau=self.get_data('CO2R_Au_kanan_normalized_RHE.csv')

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
        self.DATA_snhu=self.get_data('CO2R_Sn_Hu.csv')
        self.DATA_miyoung=self.get_data('CO2R_Miyoung.csv')

        self.DATA_bell=self.get_data('CO2R_Ag_Ezra.csv')

        self.DATA_jihun=self.get_data('CO2R_Jihun.csv')

        os.chdir('/'.join(os.path.realpath(__file__).split('/')[:-1])+'/../data/NOR')
        self.DATA_norchoi=self.get_data('NOR_choi.csv')


#        self.DATA_lr2=self.get_data('COR_lei_high_surface_x380_rhe.csv')
        os.chdir(root)
        


    def get_data(self,filename):
        o = open('CSV/'+filename,'rU')
        DATA = Object()
        DATA.data = np.array(list(csv.reader(o)))
        DATA.label=[str(d) for d in DATA.data[0,:]] #data[0,1:4]
        DATA.data = np.delete(DATA.data,0,0) # delete header
        return DATA

    def plot_stuff(self,**kwargs):
        """"
        possible kwargs:
            list_of_data,DATA,fit,joinlines, skip={}, ax=None,voltage_mode='previous',take_log=False,\
            linestyle=None,symbol=None,convert=None,fit_tafel=False,legend=False,msize=12,color=None,lw=1,fs=14,ls=None,marker=None,\
            color_mode='species'):
        """
        #if voltage_mode='previous' take previous point before list_of_data as voltages, if 'first', take first column for all data
        #if fit=0 Tafel, fit=1 poly, fit=2 poly no points
        #fit_tafel: New tafel fitting based on error estimation due to curvature
        kwargs['linestyle']='-'
        kwargs['symbol']='o'
        def assign(key):
            if key in kwargs:
                return kwargs[key]
            else:
                return None
        list_of_data=assign('list_of_data')
        DATA=assign('DATA')
        fit=assign('fit')
        joinlines=assign('joinlines')
        skip=assign('skip')
        ax=assign('ax')
        voltage_mode=assign('voltage_mode')
        take_log=assign('take_log')
        plot_mode=assign('plot_mode')
        linestyle=assign('linestyle')
        symbol=assign('symbol')
        convert=assign('convert')
        fit_tafel=assign('fit_tafel')
        legend=assign('legend')
        msize=assign('msize')
        color=assign('color')
        lw=assign('lw')
        fs=assign('fs')
        ls=assign('ls')
        marker=assign('marker')
        print('cmarker',marker)
        color_mode=assign('color_mode')

        n=0
        if ax is None:
            ax=plt
        #if take_log:
        #    func=ax.semilogy
        #else:
        #    func=ax.plot
        func=ax.plot #ax.semilogy
        markers=cycle(['o','*','x','d'])
        filled=cycle(['filled','none'])
        symbols=self.symbols
        for j in list_of_data[:-1]: #[1,3,5,7,9,11,13]:
            #elif symbol is not None:
            #    cmarker=next(symbols)
            voltage=[]
            current=[]
            for l in range(0,len(DATA.data[:,j-1])):
                try:
                    if plot_mode=='partial_current':
                        current.append(float(DATA.data[l,j]))
                    elif plot_mode=='faradaic_efficiency':
                        current.append(float(DATA.data[l,j])/float(DATA.data[l,list_of_data[-1]]))
                    if voltage_mode=='previous':
                        voltage.append(float(DATA.data[l,j-1]))
                    elif voltage_mode=='first':
                        voltage.append(float(DATA.data[l,0]))
                    else:
                        print('Invalid voltage_mode')
                except ValueError:
                    pass
#            if any(current)<0:
#                current=[-c for c in current]
            if not take_log and not plot_mode=='faradaic':
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
                print('getting',DATA.label[j])
                self.color[n]=self.get_color(DATA.label[j],mode=color_mode)
            else:
                self.color[n]=color
            val=re.findall('pH[ ]*=[ ]*([0-9.]+)',DATA.label[j])
            if len(val)>0:
                pH=float(val[0])
            else:
                pH=7.0
                print('No pH found to convert scales!! Check CSV data labels if pH is included')
                sys.exit()
            print('marker',self.marker)
            cmarker=marker
            if cmarker is None:
                cmarker=next(markers)
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
                mfc=next(filled)
                if cmarker=='*':
                    func(voltage, current, color=self.color[n], linestyle=linestyle_2, marker = cmarker, markersize=12,label=DATA.label[j],lw=lw,ls=ls,zorder=1e6) #,mfc='none')# facecolors=next(filled))  #linestyle = ':',
                else:
                    func(voltage, current, color=self.color[n], linestyle=linestyle_2, marker = cmarker, markersize=msize,label=DATA.label[j],lw=lw,ls=ls,zorder=1e6)
#                plt.show()
#                sys.exit()
            if fit_tafel:
                linestyle=':'
            if fit>0:
                if len(voltage)>1:
                    species=DATA.label[j].split(',')[0].strip()
                    skip_val=0
                    for sk in skip:
                        if sk in species:
                            skip_val=skip[sk]
                    if fit==1:
                        self.plot_tafel(ax,voltage,current,self.color[n],skip_val,linestyle=linestyle,take_log=take_log,lw=lw,fs=fs,ls=ls)
                    if fit>1:
                        V=voltage
                        J=current
                        linecol=self.color[n]
                        if 'H$_2' in DATA.label[j]:
                            order=3
                        else:
                            order=2
                        a=np.array([V,J]).T
                        a=a[a[:,0].argsort()]
                        V=a[:,0]
                        J=a[:,1]
                        J=np.log10(J)
                        spl=UnivariateSpline(V,J,k=3,s=1)
                        p=lambda x: 10**spl(x)
                        #z=np.polyfit(V,J,order)
                        #p=np.poly1d(z)
                        xfine=np.linspace(min(V),max(V),1000)
                        if linestyle=="None":
                            linestyle='-'
                        if fit>2:
                            func(xfine,p(xfine),linestyle=linestyle,color=linecol,label=DATA.label[j],lw=lw,ls=ls) #,zorder=1e7)
                        else:
                            func(xfine,p(xfine),color=linecol,lw=lw,ls=ls)
                            
            n=n+1   
        #plt.ylabel('log$j$ [mA/cm$^2_{\mathrm{real}}$] ')
        #plt.ylim((-3,1))
        #plt.xlim((-1.4,-1.05))
        if legend:
            leg = ax.legend(loc='lower left',
                  ncol=1, fontsize=6, numpoints=1) #, fancybox=True, shadow=False, loc='upper center', bbox_to_anchor=(0.5, 1.05))
            for line,text in zip(leg.get_lines(), leg.get_texts()):
                text.set_color(line.get_color())
    
    def plot_tafel(self,ax,V,J,linecol,skip=0,linestyle='-',take_log=True,lw=1,fs=11,ls=None):
        itafel+=5.
        J=np.log10(J)
        data=[V,J]
        data=list(map(list,list(zip(*data))))
        data=np.array(data)
        data.view(data.dtype.str+','+data.dtype.str).sort(order=['f0'], axis=0)
        V=data[:,0]
        J=data[:,1]
        if skip:
            if type(skip)==list:
                Vs=[V[skip[0]:skip[1]],V[skip[1]:]]
                Js=[J[skip[0]:skip[1]],V[skip[1]:]]
            else:
                Vs=[V[skip:]]
                Js=[J[skip:]]
        else:
            Vs=[V]
            Js=[J]
        for cV,cJ in zip(Vs,Js):
            if len(cV)==0 or len(cJ)==0:
                print('WARNING: skipping all data in Tafel slope plot')
                continue
            #remove all NaN
            A=np.array([[vv,jj] for vv,jj in zip(cV,cJ) if (not math.isnan(vv) and not math.isnan(jj))])
            if len(A)>0:
                V,J=A.T
            else:
                continue
            #fit data
            p=np.polyfit(V,J,1)
            z=np.poly1d(p)
            #tafel slope
            tafel_slope=-1000./p[0]
            minV=min(V)
            maxV=max(V)
            dV=maxV-minV
            Vf=np.linspace(minV-dV/4.,maxV+dV/3.,1000)
            ax.semilogy(Vf,10**z(Vf),color=linecol, linestyle=linestyle,lw=lw)
            an=ax.annotate(str(int(tafel_slope))+'mV/dec',xy=(V[0]+0.2+itafel,10**(z(V[0])-1)),color=linecol, fontsize=fs)
            an.draggable()
        return ax
    
    def plot_data(self,ax=None,reference=['all'],species=['all'],pH=['all'],ci_bic=['all'],scale='RHE',reaction='all',\
            system=['all'],coloring='species',fit_tafel=False,only_points=False,\
            take_log=True,marker=None,legend=False,msize=None,color=None,lw=1,ls=None,\
            color_mode='species',plot_mode='partial_current'):
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
        color_mode: either color the data by species ('species') or research group ('group') or metal surface ('surface')
        plot_mode: either partial_current (default) or faradaic_efficiency
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
            print('Reference, species and pH must be lists!')
            sys.exit()
    
        #work on the species list, expand C1 and C2 in the respective species:
        #C1 does by definition here NOT include HCOO (formate,formic acid)!!
        c1_list=['CH$_4$','MeOH','CH$_3$OH'] #,'HCOO','MeOH','CH$_3$OH'] #,'C$_1$']
        c2_list=['C$_2$H$_4$','CH$_3$COO','EtOH','CH$_2$CH$_2$','C$_2$H$_6$','MeCHO','GlycAld','AcetAld','EtGlyci']
        c2p_list=c2_list+['n-PrOH','PrOH','AllylAlc','Acetone']
        more_prods=['CO','H$_2$','HCOO','HCOO-','Formate']

        nor_prods=['NH$_2$OH','N$_2$O','NH$_3$']
        #contains all products
        all_prods=c1_list+c2_list+c2p_list+more_prods+nor_prods

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
    
        print(('Plotting species = {}'.format(species)))
    
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
                    print('No pH given in data file ',DATA.label)
                    print('This is the line that was searched:',line)
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
                    found=any([system==l.strip() for l in line.split(',')]) #system in [l.strip() for l in line.split(',')]#re.findall(system,line)
                    if found: #len(found)>0:
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
            if plot_mode=='faradaic_efficiency':
                total_count=0
                #append total current for plotting
                for i,h in enumerate(DATA.label):
                    if h.startswith("TOTAL") and anyph(pH,h) and reac(h) and system(h) and i not in heads:
                        total_count+=1
                        if total_count>1:
                            print('Only one system is allowed per csv file, so that only one total current column is contained, so please create a separate csv file for each system')
                            sys.exit()
                        heads+=[i]
            print('THE HEADS',heads)
            return heads
    
        linestyles=cycle(['-','-',':','-.'])
        self.symbols=cycle(['o','d','s','D','x','p','*','v','h','1','2'])
        symbols=self.symbols
    
        if fit_tafel:
            fit=1
        else:
            fit=2

        if only_points:
            fit=0
        
    
    #    data_labels=[','.join(a.label[1].split(',')[1:]) for a in [self.DATA_ss_01,self.DATA_ss_005,self.DATA_ss_02,self.DATA_ss_005_FE,self.DATA_ss_01_FE,self.DATA_ss_02_FE,self.DATA_jr,self.DATA_kr,self.DATA_lr,self.DATA_h2s,self.DATA_hr,self.DATA_hs,self.DATA_]]
        data_labels=[','.join(a.label[1].split(',')[1:]) for a in [self.DATA_ss_01,self.DATA_ss_005,self.DATA_ss_02,self.DATA_ss_005_FE,self.DATA_ss_01_FE,self.DATA_ss_02_FE,self.DATA_jr,self.DATA_jr_NF,self.DATA_jr_NC,self.DATA_kr,self.DATA_kau,self.DATA_lr,self.DATA_lr_NF,self.DATA_lr_HER,self.DATA_h2s,self.DATA_hau,self.DATA_hr,self.DATA_hs,self.DATA_wr_NW,self.DATA_wr_NW_all,self.DATA_wu,self.DATA_du,self.DATA_du2,self.DATA_cas,self.DATA_cas2,self.DATA_norchoi,self.DATA_snhu,self.DATA_miyoung,self.DATA_bell,self.DATA_jihun]]

        #basic settings applied to each experimental data set for plotting
        kwargs_base={'ax':ax,'take_log':take_log,'fit_tafel':fit_tafel,'legend':legend,'msize':msize,'lw':lw,'ls':ls,'marker':marker,'joinlines':0,'color_mode':color_mode,'fit':fit,'plot_mode':plot_mode}

        for cpH in pH:
            skip_dict={}
            for label in data_labels:
                #initialize with 0
                skip_dict[label]={}
                for prod in all_prods:
                    skip_dict[label][prod]=0

                    if 'pc-Au' in label and 'pH = 3.0' in label:
                        skip_dict[label]['CO']=5
                    elif 'pc-Au' in label and 'pH = 6.8' in label:
                        skip_dict[label]['CO']=6
                if 'Dunwell' in label:
                    if 'pc-Au' in label and 'pH = 7.3' in label:
                        skip_dict[label]['CO']=[0,7]
                elif 'Kanan' in label:
                    skip_dict[label]['HCOO']=9
                    skip_dict[label]['CO']=10
                    skip_dict[label]['H$_2$']=9
                    skip_dict[label]['C$_2$H$_4$']=2
                    skip_dict[label]['C$_2$H$_6$']=3
                elif 'Hori' in label:
                    skip_dict[label]['CO']=17
                elif 'Strasser' in label:
                    skip_dict[label]['CO']=8
            for bic in ci_bic:
                if isref('strasser'):
                    if bic in ['0.1','all']:
                        DATA=self.DATA_ss_01
                        name=','.join(DATA.label[1].split(',')[1:])
                        skip=skip_dict[name]
                        spp=s2i(species,DATA,pH=cpH)
                        if len(spp)>0:
                            linestyle=next(linestyles)
                            if color_mode=='group':
                                color=next(colors)
                            if marker is None:
                                symbol=next(symbols)
                            else:
                                symbol=marker
                            if scale=='RHE':
                                self.plot_stuff(spp,DATA,fit,0,skip,ax=ax,linestyle=linestyle,symbol=symbol,\
                                    take_log=take_log,convert='SHE_TO_RHE',fit_tafel=fit_tafel,legend=legend,\
                                    msize=msize,color=color,lw=lw,ls=ls,marker=marker)
                            elif scale=='SHE':
                                self.plot_stuff(spp,DATA,fit,0,skip,ax=ax,linestyle=linestyle,symbol=symbol,\
                                    take_log=take_log,fit_tafel=fit_tafel,legend=legend,\
                                    msize=msize,color=color,lw=lw,ls=ls,marker=marker)
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
                                    msize=msize,color=color,lw=lw,ls=ls,marker=marker)
                            elif scale=='SHE':
                                self.plot_stuff(spp,DATA,fit,0,skip,ax=ax,linestyle=linestyle,symbol=symbol,\
                                    take_log=take_log,fit_tafel=fit_tafel,legend=legend,msize=msize,color=color,lw=lw,ls=ls,marker=marker)
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
                            linestyle=next(linestyles)
                            if marker is None:
                                symbol=next(symbols)
                            else:
                                symbol=marker
                            if scale=='RHE':
                                self.plot_stuff(spp,DATA,fit,0,skip,ax=ax,linestyle=linestyle,symbol=symbol,\
                                    take_log=take_log,convert='SHE_TO_RHE',fit_tafel=fit_tafel,legend=legend,\
                                    msize=msize,color=color,lw=lw,ls=ls,marker=marker)
                            elif scale=='SHE':
                                self.plot_stuff(spp,DATA,fit,0,skip,ax=ax,linestyle=linestyle,symbol=symbol,\
                                    take_log=take_log,fit_tafel=fit_tafel,legend=legend,msize=msize,color=color,lw=lw,ls=ls,marker=marker)
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
                            kwargs=kwargs_base.copy()
                            if scale=='RHE':
                                kwargs['convert']=None #'SHE_TO_RHE'
                            elif scale=='SHE':
                                kwargs['convert']='RHE_TO_SHE'
                            kwargs['list_of_data']=spp
                            kwargs['voltage_mode']='previous'
                            kwargs['DATA']=DATA
                            kwargs['skip']=skip
                            kwargs['take_log']=True
                            self.plot_stuff(**kwargs)


                    if isref('jaramillo'):
                        DATA=self.DATA_jr
                        name=','.join(DATA.label[1].split(',')[1:])
                        skip=skip_dict[name]
                        spp=s2i(species,DATA,pH=cpH)
                        if len(spp)>0:
                            kwargs=kwargs_base.copy()
                            if scale=='RHE':
                                kwargs['convert']=None #'SHE_TO_RHE'
                            elif scale=='SHE':
                                kwargs['convert']='RHE_TO_SHE'
                            kwargs['list_of_data']=spp
                            kwargs['voltage_mode']='first'
                            kwargs['DATA']=DATA
                            kwargs['skip']=skip
                            kwargs['take_log']=True
                            self.plot_stuff(**kwargs)
                        DATA=self.DATA_lr
                        name=','.join(DATA.label[1].split(',')[1:])
                        skip=skip_dict[name]
                        spp=s2i(species,DATA,pH=cpH)
                        if len(spp)>0:
                            kwargs=kwargs_base.copy()
                            if scale=='RHE':
                                kwargs['convert']=None #'SHE_TO_RHE'
                            elif scale=='SHE':
                                kwargs['convert']='RHE_TO_SHE'
                            kwargs['list_of_data']=spp
                            kwargs['voltage_mode']='first'
                            kwargs['DATA']=DATA
                            kwargs['skip']=skip
                            kwargs['take_log']=True
                            self.plot_stuff(**kwargs)
                        DATA=self.DATA_lr_NF
                        name=','.join(DATA.label[1].split(',')[1:])
                        skip=skip_dict[name]
                        spp=s2i(species,DATA,pH=cpH)
                        if len(spp)>0:
                            kwargs=kwargs_base.copy()
                            if scale=='RHE':
                                kwargs['convert']=None #'SHE_TO_RHE'
                            elif scale=='SHE':
                                kwargs['convert']='RHE_TO_SHE'
                            kwargs['list_of_data']=spp
                            kwargs['voltage_mode']='first'
                            kwargs['DATA']=DATA
                            kwargs['skip']=skip
                            kwargs['take_log']=True
                            self.plot_stuff(**kwargs)
                        DATA=self.DATA_jr_NF
                        name=','.join(DATA.label[1].split(',')[1:])
                        skip=skip_dict[name]
                        spp=s2i(species,DATA,pH=cpH)
                        if len(spp)>0:
                            kwargs=kwargs_base.copy()
                            if scale=='RHE':
                                kwargs['convert']=None #'SHE_TO_RHE'
                            elif scale=='SHE':
                                kwargs['convert']='RHE_TO_SHE'
                            kwargs['list_of_data']=spp
                            kwargs['voltage_mode']='first'
                            kwargs['DATA']=DATA
                            kwargs['skip']=skip
                            kwargs['take_log']=True
                            self.plot_stuff(**kwargs)
                        DATA=self.DATA_jr_NC
                        name=','.join(DATA.label[1].split(',')[1:])
                        skip=skip_dict[name]
                        spp=s2i(species,DATA,pH=cpH)
                        if len(spp)>0:
                            kwargs=kwargs_base.copy()
                            if scale=='RHE':
                                kwargs['convert']=None #'SHE_TO_RHE'
                            elif scale=='SHE':
                                kwargs['convert']='RHE_TO_SHE'
                            kwargs['list_of_data']=spp
                            kwargs['voltage_mode']='first'
                            kwargs['DATA']=DATA
                            kwargs['skip']=skip
                            kwargs['take_log']=True
                            self.plot_stuff(**kwargs)
                        DATA=self.DATA_cas
                        name=','.join(DATA.label[1].split(',')[1:])
                        skip=skip_dict[name]
                        spp=s2i(species,DATA,pH=cpH)
                        if len(spp)>0:
                            #linestyle=next(linestyles)
                            #symbol=next(symbols)
                            #if scale=='SHE':
                            #    self.plot_stuff(spp,DATA,fit,0,skip,linestyle=linestyle,symbol=symbol,\
                            #        voltage_mode='previous',take_log=True,fit_tafel=fit_tafel,legend=legend,ax=ax,msize=msize,color=color,lw=lw,ls=ls,marker=marker)
                            #elif scale=='RHE':
                            #    self.plot_stuff(spp,DATA,fit,0,skip,linestyle=linestyle,symbol=symbol,\
                            #        voltage_mode='previous',take_log=True,convert='SHE_TO_RHE',fit_tafel=fit_tafel,legend=legend,\
                            #        msize=msize,ax=ax,color=color,lw=lw,ls=ls,marker=marker)
                            kwargs=kwargs_base.copy()
                            if scale=='RHE':
                                kwargs['convert']='SHE_TO_RHE'
                            elif scale=='SHE':
                                kwargs['convert']=None
                            kwargs['list_of_data']=spp
                            kwargs['voltage_mode']='previous'
                            kwargs['DATA']=DATA
                            kwargs['skip']=skip
                            self.plot_stuff(**kwargs)
                        #DATA=self.DATA_cas2
                        #name=','.join(DATA.label[1].split(',')[1:])
                        #skip=skip_dict[name]
                        #spp=s2i(species,DATA,pH=cpH)
                        #if len(spp)>0:
                        #    linestyle=next(linestyles)
                        #    symbol=next(symbols)
                        #    if scale=='SHE':
                        #        self.plot_stuff(spp,DATA,fit,0,skip,linestyle=linestyle,symbol=symbol,\
                        #            voltage_mode='previous',take_log=True,fit_tafel=fit_tafel,legend=legend,ax=ax,msize=msize,color=color,lw=lw)
                        #    elif scale=='RHE':
                        #        self.plot_stuff(spp,DATA,fit,0,skip,linestyle=linestyle,symbol=symbol,\
                        #            voltage_mode='previous',take_log=True,convert='SHE_TO_RHE',fit_tafel=fit_tafel,legend=legend,\
                        #            msize=msize,ax=ax,color=color,lw=lw)
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
                        DATA=self.DATA_wr_NW_all
                        name=','.join(DATA.label[1].split(',')[1:])
                        skip=skip_dict[name]
                        spp=s2i(species,DATA,pH=cpH)
                        if len(spp)>0:
                            kwargs=kwargs_base.copy()
                            if scale=='RHE':
                                kwargs['convert']=None #'SHE_TO_RHE'
                            elif scale=='SHE':
                                kwargs['convert']='RHE_TO_SHE'
                            kwargs['list_of_data']=spp
                            kwargs['voltage_mode']='previous'
                            kwargs['DATA']=DATA
                            kwargs['skip']=skip
                            kwargs['take_log']=True
                            self.plot_stuff(**kwargs)
                    if isref('hori'):
                        DATA=self.DATA_h2s
                        name=','.join(DATA.label[1].split(',')[1:])
                        skip=skip_dict[name]
                        spp=s2i(species,DATA,pH=cpH)
                        if len(spp)>0:
                            kwargs=kwargs_base.copy()
                            if scale=='RHE':
                                kwargs['convert']='SHE_TO_RHE'
                            elif scale=='SHE':
                                kwargs['convert']=None
                            kwargs['list_of_data']=spp
                            kwargs['voltage_mode']='first'
                            kwargs['DATA']=DATA
                            kwargs['skip']=skip
                            kwargs['take_log']=True
                            self.plot_stuff(**kwargs)
                        DATA=self.DATA_hau
                        name=','.join(DATA.label[1].split(',')[1:])
                        skip=skip_dict[name]
                        spp=s2i(species,DATA,pH=cpH)
                        if len(spp)>0:
                            kwargs=kwargs_base.copy()
                            if scale=='RHE':
                                kwargs['convert']='SHE_TO_RHE'
                            elif scale=='SHE':
                                kwargs['convert']=None
                            kwargs['list_of_data']=spp
                            kwargs['voltage_mode']='previous'
                            kwargs['DATA']=DATA
                            kwargs['skip']=skip
                            kwargs['take_log']=False
                            self.plot_stuff(**kwargs)
    
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
                            self.plot_stuff(spp,DATA,fit,0,skip,ax=ax,linestyle=linestyle,symbol=symbol,fit_tafel=fit_tafel,legend=legend,msize=msize,color=color,lw=lw,ls=ls,marker=marker)
    
                    if isref('wuttig'):
                        DATA=self.DATA_wu
                        name=','.join(DATA.label[1].split(',')[1:])
                        skip=skip_dict[name]
                        spp=s2i(species,DATA,pH=cpH)
                        if len(spp)>0:
                            kwargs=kwargs_base.copy()
                            if scale=='RHE':
                                kwargs['convert']='SHE_TO_RHE'
                            elif scale=='SHE':
                                kwargs['convert']=None
                            kwargs['list_of_data']=spp
                            kwargs['voltage_mode']='first'
                            kwargs['DATA']=DATA
                            kwargs['skip']=skip
                            kwargs['take_log']=False
                            self.plot_stuff(**kwargs)

    
                    if isref('choi'):
                        DATA=self.DATA_norchoi
                        name=','.join(DATA.label[1].split(',')[1:])
                        skip=skip_dict[name]
                        spp=s2i(species,DATA,pH=cpH)
                        if len(spp)>0:
                            kwargs=kwargs_base.copy()
                            if scale=='RHE':
                                kwargs['convert']=None
                            elif scale=='SHE':
                                kwargs['convert']='RHE_TO_SHE'
                            kwargs['list_of_data']=spp
                            kwargs['voltage_mode']='first'
                            kwargs['DATA']=DATA
                            kwargs['skip']=skip
                            kwargs['take_log']=True
                            self.plot_stuff(**kwargs)

                    if isref('hu'):
                        DATA=self.DATA_snhu
                        name=','.join(DATA.label[1].split(',')[1:])
                        skip=skip_dict[name]
                        spp=s2i(species,DATA,pH=cpH)
                        if len(spp)>0:
                            kwargs=kwargs_base.copy()
                            if scale=='RHE':
                                kwargs['convert']=None
                            elif scale=='SHE':
                                kwargs['convert']='RHE_TO_SHE'
                            kwargs['list_of_data']=spp
                            kwargs['voltage_mode']='first'
                            kwargs['DATA']=DATA
                            kwargs['skip']=skip
                            kwargs['take_log']=True
                            self.plot_stuff(**kwargs)

                    if isref('miyoung'):
                        DATA=self.DATA_miyoung
                        name=','.join(DATA.label[1].split(',')[1:])
                        skip=skip_dict[name]
                        spp=s2i(species,DATA,pH=cpH)
                        if len(spp)>0:
                            kwargs=kwargs_base.copy()
                            if scale=='RHE':
                                kwargs['convert']='SHE_TO_RHE'
                            elif scale=='SHE':
                                kwargs['convert']=None
                            kwargs['list_of_data']=spp
                            kwargs['voltage_mode']='previous'
                            kwargs['DATA']=DATA
                            kwargs['skip']=skip
                            kwargs['take_log']=False
                            self.plot_stuff(**kwargs)


                    if isref('jihun'):
                        DATA=self.DATA_jihun
                        name=','.join(DATA.label[1].split(',')[1:])
                        skip=skip_dict[name]
                        spp=s2i(species,DATA,pH=cpH)
                        if len(spp)>0:
                            kwargs=kwargs_base.copy()
                            if scale=='RHE':
                                kwargs['convert']=None #'SHE_TO_RHE'
                            elif scale=='SHE':
                                kwargs['convert']='RHE_TO_SHE'
                            kwargs['list_of_data']=spp
                            kwargs['voltage_mode']='previous'
                            kwargs['DATA']=DATA
                            kwargs['skip']=skip
                            kwargs['take_log']=True
                            self.plot_stuff(**kwargs)


                    if isref('bell'):
                        DATA=self.DATA_bell
                        name=','.join(DATA.label[1].split(',')[1:])
                        skip=skip_dict[name]
                        spp=s2i(species,DATA,pH=cpH)
                        if len(spp)>0:
                            kwargs=kwargs_base.copy()
                            if scale=='RHE':
                                kwargs['convert']='SHE_TO_RHE'
                            elif scale=='SHE':
                                kwargs['convert']=None
                            kwargs['list_of_data']=spp
                            kwargs['voltage_mode']='previous'
                            kwargs['DATA']=DATA
                            kwargs['skip']=skip
                            kwargs['take_log']=False
                            self.plot_stuff(**kwargs)

                    if isref('dunwell'):
                        DATA=self.DATA_du
                        name=','.join(DATA.label[1].split(',')[1:])
                        skip=skip_dict[name]
                        spp=s2i(species,DATA,pH=cpH)
                        if len(spp)>0:
                            kwargs=kwargs_base.copy()
                            if scale=='RHE':
                                kwargs['convert']=None
                            elif scale=='SHE':
                                kwargs['convert']='RHE_TO_SHE'
                            kwargs['list_of_data']=spp
                            kwargs['voltage_mode']='first'
                            kwargs['DATA']=DATA
                            kwargs['skip']=skip
                            kwargs['take_log']=False
                            self.plot_stuff(**kwargs)
                        DATA=self.DATA_du2
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
                                    symbol=symbol,fit_tafel=fit_tafel,legend=legend,msize=msize,color=color,lw=lw,ls=ls,marker=marker)
                            elif scale=='SHE':
                                self.plot_stuff(spp,DATA,fit,0,skip,ax=ax,voltage_mode='first',take_log=False,linestyle=linestyle,\
                                    symbol=symbol,convert='RHE_TO_SHE',fit_tafel=fit_tafel,legend=legend,msize=msize,color=color,lw=lw,ls=ls,marker=marker)

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
                                    symbol=symbol,fit_tafel=fit_tafel,legend=legend,msize=msize,color=color,lw=lw,ls=ls,marker=marker)
                            elif scale=='SHE':
                                self.plot_stuff(spp,DATA,fit,0,skip,ax=ax,voltage_mode='first',take_log=take_log,linestyle=linestyle,\
                                    symbol=symbol,convert='RHE_TO_SHE',fit_tafel=fit_tafel,legend=legend,msize=msize,color=color,lw=lw,ls=ls,marker=marker)
                        DATA=self.DATA_kau
                        name=','.join(DATA.label[1].split(',')[1:])
                        skip=skip_dict[name]
                        spp=s2i(species,DATA,pH=cpH)
                        if len(spp)>0:
                            kwargs=kwargs_base.copy()
                            if scale=='RHE':
                                kwargs['convert']=None
                            elif scale=='SHE':
                                kwargs['convert']='RHE_TO_SHE'
                            kwargs['DATA']=DATA
                            kwargs['list_of_data']=spp
                            kwargs['voltage_mode']='previous'
                            kwargs['skip']=skip
                            self.plot_stuff(**kwargs)

    def get_color(self,species,mode='species'):
        """
        mode='species':
            define color associated with a particular species
        mode='group':
            define color for each research group dataset
        mode='surface':
            define color for each investigated surface/catalyst/electrode
        """
        if mode=='species':
            species=species.strip()
            color=None #'k'
            if any([species.startswith(a) for a in ['n-PrOH']]):
                color='lightblue'
            elif any([species.startswith(a) for a in ['CH$_3$CH$_2$OH','CH3CH2OH','C$_2$H$_4$','CH$_3$COO','etol','C2','EtOH','C$_{2+}$','C$_2$','NH$_2$OH','NH2OH','C2+-sum','C2-sum']]):
                color='r'
            elif any([species.startswith(a) for a in ['CH$_4$','C$_1$','CH4','C1','N$_2$O','N2O']]):
                color='orange'
            elif any([species.startswith(a) for a in ['HCOO','HCOO-','HCOOH','Formate',"NH$_3$",'NH3']]):
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
        elif mode=='surface':
            species=species.strip()
            if 'SnO/C' in species:
                color='C0'
            elif 'SnO2/C' in species:
                color='C1'
            elif 'SnO-I/C' in species:
                color='C2'
            elif 'Sn/C' in species:
                color='C3'
            elif 'pc-Au' in species:
                color='C4'
            elif 'pc-Ag' in species:
                color='C5'
            elif 'pc-Cu' in species:
                color='C6'
            elif 'Ag' in species:
                color='C7'
            else:
                print('No color assigned for this surface, add it in experimental.get_color')
                sys.exit
        elif mode=='group':
            species=species.strip()
            if 'Dunwell' in species:
                color='C0'
            elif 'Wuttig' in species:
                color='C1'
            elif 'Jaramillo' in species:
                color='C2'
            elif 'Hori' in species:
                color='C3'
            elif 'Kanan' in species:
                color='C4'
            elif 'Hu' in species:
                color='C5'
            elif 'Bell' in species:
                color='C6'
            elif 'Jihun' in species:
                color='C7'
            else:
                print('No color assigned for this research group, add it in experimental.get_color')
                sys.exit
        else:
            print('No color_mode like {}'.format(mode))
            sys.exit()
        print(species,mode)
        return color
                #color[n]='dark red'
