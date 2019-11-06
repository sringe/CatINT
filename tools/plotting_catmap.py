import re
from tabulate import tabulate
import os
import matplotlib.pyplot as plt
from itertools import cycle
import sys
#sys.path.insert(0,'/scratch/users/sringe/transport/CatINT')
sys.path.insert(0,'/scratch/users/sringe/transport/catmap')
sys.path.insert(0,'/scratch/users/sringe/transport/catint2')
import numpy as np
from glob import glob
from catint.experimental import EXPDATA
from matplotlib import rc
from catint.io import reglob
import argparse
from catmap import ReactionModel

parser = argparse.ArgumentParser(description='Process some arguments.')

parser.add_argument('--file',help='file names to be evaluated',nargs='+')
parser.add_argument('--elemrates',help='If true show the full rate decomposition',action="store_true")
parser.add_argument('--products',help='Plot only data of selected product',nargs='+')
parser.add_argument('--coverages',help='Plot only data of selected product',nargs='+')
parser.add_argument('--pH',help='System pH. If not given, this will be tried to be read from the catmap input file, otherwise set to 7')
parser.add_argument('--expdata',help='plot also experimental data',action="store_true")
parser.add_argument('--systems',help='List of systems, for which exp. data is plotted: all, pc-Cu (default), nf-Cu, nc-Cu',nargs='+')
parser.add_argument('--scale',help='RHE (default) or SHE')
parser.add_argument('--color',help='pH or species (default)')
parser.add_argument('--title',help='Title of figure')
parser.add_argument('--name',help='File name of figure')
parser.add_argument('--fit_tafel',help='If tafel slope should be fitted to experimental data',action="store_true")
parser.add_argument('--exp_colormode',help='colormode, either color according to species or dataset')
parser.add_argument('--exp_add_pH',help='parse list of additional pH values for which experimental data should be plotted',nargs='+')
parser.add_argument('--ratecontrol',help='plot also degree of rate-control',action="store_true")

args = parser.parse_args()

msize=8
lw=2

skip=1

if args.scale is None:
    args.scale='RHE'

if args.color is None:
    args.color='species'

if args.title is None:
    args.title='catmap_fig'

if args.name is None:
    args.name='catmap_fig'

if args.exp_add_pH is None:
    args.exp_add_pH=[]

if args.exp_colormode is None:
    args.exp_colormode='species'


#if args.products is None:
#    args.products='all'

color_pH={
    '3.0':'C5',
    '6.0':'C1',
    '6.5':'C1',
    '6.8':'C2',
    '7.0':'C4',
    '7.15':'C3',
    '7.2':'C3',
    '8.0':'C2',
    '9.0':'C2',
    '10.0':'C2',
    '12.2':'C3',
    '13.0':'C3'}

show_legend=True

if args.systems is None:
    systems=['pc-Cu']
else:
    systems=args.systems

#first read all variables from one of the mkm files

#rc('text', usetex=False)

exp=EXPDATA()

j_log_plot=True

fig=plt.figure(figsize=(5, 7))
#fig=plt.figure(figsize=(3, 5))

colorlist={}
colorlist['CO']='r'
colorlist['H2']='b'
colorlist['CO2']='k'
colorlist['CH4']='olive'
colorlist['CH3CH2OH']='orange'
colorlist['H']='lightblue'
colorlist['CHOH']='y'
colorlist['OCCOH']='darkred'
colorlist['CHO']='0.75'
colorlist['OCCO']='g'
colorlist['COOH']='darkblue'
symbols={
        'CO':'x',
        'H2':'o',
        'COOH':'d',
        'HCOOH':'3',
        'CO2':'D',
        'CH4':'p',
        'CH3CH2OH':'1',
        'CHO':'2'
        }



def plot_leis_new_data(ax):
    #data=np.loadtxt('her_pcCu_lei.csv')
    i=0
    for f in ['her_pcCu_lei.csv','her_pcCu_lei_NF_norm.csv','her_pcCu_lei_NF.csv']:
        i+=1
        data=np.loadtxt(f)
        ax.semilogy(data[:,0],10**(data[:,1]),str(i+1),color='k')
    return ax

def read_data(files,header=False,dtype='species'):
    #dtype=species or elem (elementary reaction)
    pdata={}
    iarg=-1
    labels=[]
    for arg2 in files:
        iarg+=1
        header_txt=None
        #product selection
        if (dtype == 'species' and args.products is not None) or (dtype=='cov' and args.coverages is not None):
            ads=arg2.split('/')[-1].split('_')[1].split('.')[0]
            if (dtype=='species' and ads not in args.products) or (dtype=='cov' and ads not in args.coverages):
                continue
        if dtype=='elem' and args.products is not None:
            #elementary reaction index of current datafile
            inx=int(arg2.split('/')[-1].split('_')[2].split('.')[0])
            #get all product names in which this step plays a role
            names=[]
            for sp in model.rxn_mechanisms:
                indices=model.rxn_mechanisms[sp]
                if inx in indices:
                    names.append(sp.split('_')[0])
            #check if the desired product is any of these 
            if not any([n in args.products for n in names]):
                continue
            else:
                tmp=model.rxn_expressions[inx-1].split('<->')
                tmp_txt=''
                for ii in range(2):
                    if ii==0:
                        ik=0
                    else:
                        tmp_txt+=' $>$ '
                        ik=-1
                    fs_txt=tmp[ik].split(';')[0].replace('*','*').replace('_','\_').split('+')
                    fs_txt_new=''
                    for f in fs_txt:
                        if f.strip() not in ['ele\_g','H\_g','OH\_g','*\_t','*\_s','H2O\_g','*\_dl','*\_g','2*\_t']:
                            tmp3=f.strip().split('_')
                            if tmp3[-1] in ['g','s','dl','t']:
                                tmp2='_'.join(tmp3[:-1]).strip()
                            else:
                                tmp2=f
                            fs_txt_new+=tmp2.rstrip('\\')
                    tmp_txt+=fs_txt_new
                header_txt=tmp_txt
                #if not 'CO2' in header_txt and 'COOH' in header_txt:
                #    continue
                
        if header and header_txt is None:
            header_txt=arg2.split('_')[-1].split('.')[0]
#        with open(arg2,'r') as f:
#            header_txt=[line for iline,line in enumerate(f) if iline==0 and header]
#            if len(header_txt)>0:
#                header_txt=header_txt[0]
#            else:
#                header_txt=None
#            if header:
#                skiprows=1
#            else:
#                skiprows=0
        data=np.loadtxt(arg2) #, skiprows=skiprows)
        labels.append(header_txt)
        if len(np.shape(data))>1:
            x=data[:,0]
            y=data[:,1]
        else:
            x=data[0]
            y=data[1]
        species=arg2.split('/')[-1].split('_')[1].split('.')[0]
        if args.ratecontrol and len(arg2.split('/')[-1].split('_'))>2:
            ratecontrol=True
        else:
            ratecontrol=False
        if ratecontrol:
            species2=arg2.split('/')[-1].split('_')[2].split('.')[0]
        if species not in pdata:
            if not ratecontrol:
                pdata[species]=[]
            else:
                pdata[species]={}
        if ratecontrol and species2 not in pdata[species]:
            pdata[species][species2]=[]
        if len(np.shape(data))>1:
            if not ratecontrol:
                pdata[species].append([x[0],y[0]])
            else:
                pdata[species][species2].append([x[0],y[0]])
        else:
            if not ratecontrol:
                pdata[species].append([x,y])
            else:
                pdata[species][species2].append([x,y])
    if not ratecontrol:
        for sp in pdata:
            pdata[sp].sort(key=lambda x: x[0])
            pdata[sp]=np.array(pdata[sp])
    else:
        for sp in pdata:
            for sp2 in pdata[sp]:
                pdata[sp][sp2].sort(key=lambda x: x[0])
                pdata[sp][sp2]=np.array(pdata[sp][sp2])
    return pdata, labels

colors=cycle(['C'+str(i) for i in range(10)])
color_eps={6:'C1',10:'C2',80:'C3'}
linestyles=cycle(['-','--',':','-.'])
symbols=cycle(['o','d','1','2','x'])
k=-1
all_pH=[]
symbol_pH={}

ax1=plt.subplot('211')
if args.ratecontrol:
    ax1m=ax1.twinx()
ax2=plt.subplot('212')

colors_rc={}
for arg in args.file: #sys.argv[1:]:
    k+=1
    results_folder=arg+'/catmap_output'
    #if True: # k==2 or k==0:
    #    linestyle=next(linestyles)
    if 'etp' in arg:
        linestyle='--'
    else:
        linestyle='-'
    symbol=next(symbols)
    
    
    ax1.set_title('Polarization')
    searchedfile=[s for s in glob(arg+'/catmap_input/*/*.mkm') if 'template' not in s]
    try:
        mkm_file = sorted( searchedfile, key = lambda file: os.path.getctime(file))[-1]
        print('Reading model = {}'.format(mkm_file))
        model = ReactionModel(setup_file = mkm_file)
    except NameError:
        mkm_file = sorted( searchedfile, key = lambda file: os.path.getctime(file))[-3]
        print('Reading model = {}'.format(mkm_file))
        model = ReactionModel(setup_file = mkm_file)
    #try to read pH from catmap input file:
    if args.pH is None:
        try:
            print('Reading pH from bulk_ph = {}'.format(model.bulk_ph))
            pH=model.bulk_ph
        except:
            print('Could not find bulk pH, please specify at input')
            sys.exit()
    #    pH = None
    #    for line in open(catmap_inputs,'r'):
    #        if line.lstrip().startswith('pH'):
    #            pH=float(line.split('=')[-1])
    #            print 'Found pH in catmap input file = ',pH
    #    if pH is None:
    #        print 'Found no pH in catmap input file, assuming pH = 7'
    #        pH = 7
    else:
        pH=float(args.pH)
    print('The pH is {}'.format(pH))
    all_pH.append(pH)
    symbol_pH[pH]=symbol
    if color_pH is not None:
        try:
            color=color_pH[str(pH)]
        except:
            color='k'
#    if args.scale=='SHE':
#        pH=0.

    #READ CURRENT DENSITIES
    pdata,dummy=read_data(glob(results_folder+'/*/j_*'))
    if args.ratecontrol:
        rc_files=glob(results_folder+'/*/rc_*')
        if len(rc_files)>0:
            ratecontrol=True
            pdata_rc,dummy=read_data(rc_files)
        else:
            ratecontrol=False
    else:
        ratecontrol=False
    if args.elemrates:
        pdata_elem,label_elem=read_data(glob(results_folder+'/*/jelem*'),header=True,dtype='elem')
    cdata,dummy=read_data(glob(results_folder+'/*/cov*'),dtype='cov')
    cov_data=[]
    for isp,sp in enumerate(pdata):
        if args.products is not None and sp not in args.products:
            continue
        x=pdata[sp][:,0]
        y=pdata[sp][:,1]
        if ratecontrol:
            xrc={}
            yrc={}
            for sp2 in pdata_rc[sp]:
                xrc[sp2]=pdata_rc[sp][sp2][:,0]
                yrc[sp2]=pdata_rc[sp][sp2][:,1]
        if args.color=='species':
            color=exp.get_color(sp,mode='species')
            if color is None and sp in colorlist:
                color=colorlist[sp]
            elif color is None:
                color='k'
        elif args.color=='pH':
            if color_pH is None:
                color=exp.get_color(sp)
                if color is None and sp in colorlist:
                    color=colorlist[sp]
                elif color is None:
                    color='k'
        sol=re.findall('eps([0-9.]+)',results_folder)
        if len(sol)>0:
            eps=int(sol[0])
        else:
            eps=6.0
#        color=color_eps[eps]
        #symbol=''
        #if isp==len(pdata)-1:
        #    symbol='o'
        if 'etp' in arg:
            symbol=''
        elif 'w_tp' in arg:
            symbol='o'
            msize=3
        else:
            symbol='*'
            msize=4
        if j_log_plot:
            func=ax1.semilogy
        else:
            func=ax1.plot
#        symbol='' #next(symbols)
        if args.scale=='SHE':
            pHtmp=pH
            pH=0
        if k==0:
            func(x[skip:]+0.059*pH,y[skip:],linestyle+symbol,color=color,label=sp+', CatINT',ms=msize) #,label=arg.split('/')[-1])
        else:
            func(x[skip:]+0.059*pH,y[skip:],linestyle+symbol,color=color,ms=msize)
        if ratecontrol:
            for sp2 in pdata_rc[sp]:
                if sp2 not in colors_rc:
                    colors_rc[sp2]=next(colors)
                color_rc=colors_rc[sp2]
                ax1m.plot(xrc[sp2][skip:]+0.059*pH,yrc[sp2][skip:],':',color=color_rc,lw=1.5,label=sp2)
                ax1m.annotate(sp2,xy=(xrc[sp2][len(xrc[sp2][skip:])/3],yrc[sp2][len(xrc[sp2][skip:])/3]),color=color_rc,fontsize=12)
        if args.scale=='SHE':
            pH=pHtmp
    if args.elemrates:
        #first read in the rxn mechanism:
        isp=-1
        for sp,label in zip(pdata_elem,label_elem):
            isp+=1
            x=pdata_elem[sp][:,0]
            y=pdata_elem[sp][:,1]
            color=exp.get_color(sp)
            if color is None and sp in colorlist:
                color=colorlist[sp]
            elif color is None:
                color='k'
            if args.color=='species':
                color=exp.get_color(sp)
                if color is None and sp in colorlist:
                    color=colorlist[sp]
                elif color is None:
                    color='k'
            elif args.color=='pH':
                color=next(colors) #color_pH[str(pH)]
            symbol=''
            linestyle='--'
            #if isp==len(pdata)-1:
            #    symbol='o'
            #linestyle=next(linestyles)
            symbol=''
            if j_log_plot:
                func=ax1.semilogy
            else:
                func=ax1.plot
            if args.scale=='SHE':
                pHtmp=pH
                pH=0
            func(x[skip:]+0.059*pH,y[skip:],linestyle+symbol,color=color,label=label) #,label=arg.split('/')[-1])
            if args.scale=='SHE':
                pH=pHtmp
    for isp,sp in enumerate(cdata):
        x=cdata[sp][:,0]
        y=cdata[sp][:,1]
        color=exp.get_color(sp)
        if color is None and sp in colorlist:
            color=colorlist[sp]
        elif color is None:
            color='k'
        #if sp in symbols:
        #    symbol=symbols[sp]
        #else:
        #    symbol=''
        #if k==0:
        #    symbol='x'
        #else:
        #    symbol='o'
        if 'etp' in arg:
            symbol=''
        else:
            symbol='o'
            msize=3
        #if isp==len(cdata)-1:
        #    symbol='o'
        #linestyle=next(linestyles)
        if args.scale=='SHE':
            pHtmp=pH
            pH=0
        func=ax2.semilogy
        if k==0:
            func(x[skip:]+0.059*pH,y[skip:],linestyle+symbol,color=color,label=sp,lw=lw,ms=msize)
        else:
            func(x[skip:]+0.059*pH,y[skip:],linestyle+symbol,color=color,lw=lw,ms=msize)
        cov_data.append([x[skip:]+0.059*pH,y[skip:]])
        if args.scale=='SHE':
            pH=pHtmp
    cov_data=np.array(cov_data)
    x=cov_data[0,0,:]
    y=np.sum(cov_data[:,1,:],axis=0)
    func(x,1-y,linestyle+symbol,color='0.5',ms=msize,lw=lw,label='*')
    #sys.exit()
    if show_legend:
        ax1.legend(fontsize=10)
#    data=np.loadtxt('data.txt')
#    ax1.plot(data[:,0],data[:,1],'-')
    ax2.set_title('Coverages')
#    iarg=-1
#    for arg2 in glob(results_folder+'/*/cov*'):
#        data=np.loadtxt(arg2)
#        x=data[:,0]
#        y=data[:,1]
#        species=arg2.split('/')[-1].split('_')[-1].split('.')[0]
#        if k==0:
#            ax2.semilogy(x,y,'o',color=colorlist[species],label=arg2.split('/')[-1])
#        else:
#            ax2.semilogy(x,y,'-',color=colorlist[species])
    if True: #show_legend:
        ax2.legend(fontsize=8,ncol=2)
all_prods=['H$_2$','CO','CH$_4$','EtOH'] #'C2-sum','HCOOH']
def name_to_cm(name):
    if name=='H2':
        return 'H$_2$'
    elif name in ['CH4','CH$_4$']:
        return 'CH$_4$'# 1-sum'
    elif name in ['EtOH','CH3CH2OH','CH2CH2']:
        return 'C2-sum'
    else:
        return name
if args.products is not None:
    all_prods=[name_to_cm(a) for a in args.products]
else:
    all_prods=[name_to_cm(a) for a in all_prods]

for pH in set(all_pH+[float(a) for a in args.exp_add_pH]):
    if args.expdata:
        symbol=None #symbol_pH[pH]
        if args.color=='pH':
            color=color_pH[str(pH)]
        elif args.color=='species':
            color=None
        refs=['all']
        exp.plot_data(reference=refs,ax=ax1,species=all_prods,pH=[pH],\
            system=systems,scale=args.scale,only_points=True,\
            take_log=j_log_plot,marker=symbol,legend=show_legend,msize=3,color=color,\
            fit_tafel=args.fit_tafel,color_mode=args.exp_colormode)
#        if pH == 13:
#            #pure HER
#            exp.plot_data(reference=['jaramillo-her'],ax=ax1,species=all_prods,pH=['13.0'],\
#                system=systems,scale=args.scale,only_points=True,\
#                take_log=j_log_plot,marker=symbol,legend=show_legend,msize=3,color=color)
#            #HER + COR
#            exp.plot_data(reference=['hori','jaramillo','wang'],ax=ax1,species=all_prods,pH=['13.0'],\
#                system=systems,scale=args.scale,only_points=True,\
#                take_log=j_log_plot,marker=symbol,legend=show_legend,msize=3,color=color)
#        elif pH == 6.8 or pH == 7.0 or pH == 7.2:
#            #exp.plot_data(reference=['hori','jaramillo','wang'],ax=ax1,species=all_prods,pH=['6.8','7.0','7.2'],\
#            #    system=systems,scale=args.scale,only_points=True,\
#            #    take_log=j_log_plot,marker=symbol,legend=show_legend,msize=3,color=color)
#            only_points=True
#            fit_tafel=True
#            #wuttig
#            refs=['jaramillo'] #,'dunwell'] #,'wuttig']
#            exp.plot_data(reference=refs,ax=ax1,species=all_prods,pH=[str(pH)],\
#                system=systems,scale=args.scale,only_points=only_points,\
#                take_log=j_log_plot,marker=symbol,legend=show_legend,msize=5,color=color,fit_tafel=fit_tafel)
#            fit_tafel=False
##            exp.plot_data(reference=refs,ax=ax1,species=all_prods,pH=['3.0'],\
##                system=systems,scale=args.scale,only_points=only_points,\
##                take_log=j_log_plot,marker=symbol,legend=show_legend,msize=5,color=color,fit_tafel=fit_tafel)
#        elif pH == 3.0: #6.8 or pH == 7.0:
#            #exp.plot_data(reference=['hori','jaramillo','wang'],ax=ax1,species=all_prods,pH=['6.8','7.0','7.2'],\
#            #    system=systems,scale=args.scale,only_points=True,\
#            #    take_log=j_log_plot,marker=symbol,legend=show_legend,msize=3,color=color)
#            exp.plot_data(reference=['jaramillo','wuttig'],ax=ax1,species=all_prods,pH=['3.0'],\
#                system=systems,scale=args.scale,only_points=True,\
#                take_log=j_log_plot,marker=symbol,legend=show_legend,msize=3,color=color,fit_tafel=fit_tafel)
#        else:
#            for i,p in enumerate(all_prods):
#                if p=='CH$_4$':
#                    all_prods[i]='C$_1$'
#                elif p=='EtOH':
#                    all_prods
#            exp.plot_data(reference=['hori'],ax=ax1,species=all_prods,pH=[str(pH)],\
#                system=systems,scale=args.scale,only_points=True,\
#                take_log=j_log_plot,marker=symbol,legend=show_legend,msize=3,color=color)
#            #exp.plot_data(reference=['strasser'],ax=ax1,species=all_prods,pH=[str(pH)],\
#            #    system=systems,scale=args.scale,only_points=True,\
#            #    take_log=j_log_plot,marker=symbol,legend=show_legend,msize=3,color=color)
#ax1=plot_leis_new_data(ax1)
ax1.set_ylim([1e-5,1e2])

#def name_to_title(name):
#    if name=='CHO-H2O-ele*_t':
#        return 1, 'CHO to CHOH'
#    elif name=='CH-OH-ele*_t':
#        return 2, 'CHOH to CH'
#    elif name=='H2O-CO-ele*_t':
#        return 0, 'CO to CHO'
#
#rprint=[[0,0,0],[0,0,0]]
#for a in args.title.split(';'):
#    tmp=a.split() #[aa.strip for aa in a.split()]
#    if len(tmp)<1:
#        continue
#    inx,name=name_to_title(tmp[0])
#    if tmp[1]=='energy':
#        inx2=0
#    else:
#        inx2=1
#    rprint[inx2][inx]=tmp[2]
#
#print 'rprint = ',rprint
#pp=''
#pp+=args.name.replace('_','\_')+'\n'
#for rr in rprint:
#    if rr==rprint[0]:
#        pp+='energy '
#    else:
#        pp+='beta '
#    pp+=' '.join(rr)
#    if rr !=rprint[-1]:
#        pp+='\n'
#print 'TESTING',pp
#ax1.set_title(pp,fontsize=12)
#ax1.set_title(args.title.replace('_','-').replace('*',''),fontsize=6) #[a+'\n' for a in args.title.split(';')],fontsize=6)
if args.scale=='RHE':
    ax1.set_xlabel(r'Voltage vs. RHE (V)')
    ax2.set_xlabel(r'Voltage vs. RHE (V)')
else:
    ax1.set_xlabel(r'Voltage vs. SHE (V)')
    ax2.set_xlabel(r'Voltage vs. SHE (V)')
ax1.set_ylabel('Current Density (mV/dec)')
ax2.set_ylabel('Coverage')
plt.tight_layout()
print('Saving to {}'.format(args.name+'.pdf'))
plt.savefig(args.name+'.pdf')# ,dpi=500)
    
plt.show()
