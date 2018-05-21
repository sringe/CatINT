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

args = parser.parse_args()

print args.products

#first read all variables from one of the mkm files

#rc('text', usetex=False)

exp=EXPDATA()

j_log_plot=True

fig=plt.figure(figsize=(6, 8))

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


products=['CO','CH4','CH3CH2OH','H2','HCOOH']

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
            print 'checking inx = ',inx
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
        species=arg2.split('/')[-1].split('_')[-1].split('.')[0]
        if species not in pdata:
            pdata[species]=[]
        if len(np.shape(data))>1:
            pdata[species].append([x[0],y[0]])
        else:
            pdata[species].append([x,y])
    for sp in pdata:
        pdata[sp].sort(key=lambda x: x[0])
        pdata[sp]=np.array(pdata[sp])
    return pdata, labels

colors=cycle(['C'+str(i) for i in range(10)])
linestyles=cycle(['-','--',':','-.'])
symbols=cycle(['o','d','1','2','x'])
k=-1
all_pH=[]
symbol_pH={}
for arg in args.file: #sys.argv[1:]:
    k+=1
    results_folder=arg+'/catmap_output'
    linestyle=next(linestyles)
    symbol=next(symbols)
    
    ax1=plt.subplot('211')
    ax2=plt.subplot('212')
    
    ax1.set_title('Polarization')
    searchedfile=glob(arg+'/catmap_input/*/*.mkm')
    try:
        mkm_file = sorted( searchedfile, key = lambda file: os.path.getctime(file))[-1]
        model = ReactionModel(setup_file = mkm_file)
    except NameError:
        mkm_file = sorted( searchedfile, key = lambda file: os.path.getctime(file))[-3]
        model = ReactionModel(setup_file = mkm_file)
    
    #try to read pH from catmap input file:
    if args.pH is None:
        pH=model.pH
    #    pH = None
    #    for line in open(catmap_inputs,'r'):
    #        if line.lstrip().startswith('pH'):
    #            pH=float(line.split('=')[-1])
    #            print 'Found pH in catmap input file = ',pH
    #    if pH is None:
    #        print 'Found no pH in catmap input file, assuming pH = 7'
    #        pH = 7
    else:
        pH=args.pH
    all_pH.append(pH)
    symbol_pH[pH]=symbol
    pdata,dummy=read_data(glob(results_folder+'/*/j_*'))
    if args.elemrates:
        pdata_elem,label_elem=read_data(glob(results_folder+'/*/jelem*'),header=True,dtype='elem')
    cdata,dummy=read_data(glob(results_folder+'/*/cov*'),dtype='cov')
    for isp,sp in enumerate(pdata):
        if sp not in products:
            continue
        x=pdata[sp][:,0]
        y=pdata[sp][:,1]
        color=exp.get_color(sp)
        if color is None and sp in colorlist:
            color=colorlist[sp]
        elif color is None:
            color='k'
        symbol=''
        if isp==len(pdata)-1:
            symbol='o'
        if j_log_plot:
            func=ax1.semilogy
        else:
            func=ax1.plot
        symbol='' #next(symbols)
        if k==0:
            func(x+0.059*pH,y,linestyle+symbol,color=color,label=sp) #,label=arg.split('/')[-1])
        else:
            func(x+0.059*pH,y,linestyle+symbol,color=color)
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
            color=next(colors) #str(0.2+isp*0.2)
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
            func(x+0.059*pH,y,linestyle+symbol,color=color,label=label) #,label=arg.split('/')[-1])
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
        symbol=''
        #if isp==len(cdata)-1:
        #    symbol='o'
        print 'the color',k,sp,color
        #linestyle=next(linestyles)
        if k==0:
            ax2.plot(x+0.059*pH,y,linestyle+symbol,color=color,label=sp)
        else:
            ax2.plot(x+0.059*pH,y,linestyle+symbol,color=color)
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
    ax2.legend(fontsize=8,ncol=2)
all_prods=['H$_2$','CO','CH$_4$','EtOH'] #'C2-sum','HCOOH']
def name_to_cm(name):
    if name=='H2':
        return 'H$_2$'
    elif name=='CH4':
        return 'CH$_4$'
    elif name=='CH3CH2OH':
        return 'C2-sum'
    else:
        return name
if args.products is not None:
    all_prods=[name_to_cm(a) for a in args.products]
for pH in set(all_pH):
    if args.expdata:
        symbol=symbol_pH[pH]
        if pH == 13:
            exp.plot_data(reference=['hori','jaramillo'],ax=ax1,species=all_prods,pH=['13.0'],system=['pc-Cu'],scale='RHE',only_points=True,take_log=j_log_plot,marker=symbol)
        elif pH == 6.8 or pH == 7.0:
            exp.plot_data(reference=['hori','jaramillo'],ax=ax1,species=all_prods,pH=['6.8','7.0','7.2'],system=['pc-Cu'],scale='RHE',only_points=True,take_log=j_log_plot,marker=symbol)
#ax1=plot_leis_new_data(ax1)
ax1.set_ylim([1e-5,1e2])
ax1.set_xlim([-1.6,0.1])
ax2.set_xlim([-1.6,0.1])
ax1.set_xlabel(r'Voltage vs. RHE (V)')
ax2.set_xlabel(r'Voltage vs. RHE (V)')
plt.tight_layout()
plt.show()
