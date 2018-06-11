import re
from itertools import cycle
import os
from glob import glob
import numpy as np
import matplotlib.pyplot as plt
import sys
from catmap import ReactionModel
import argparse

parser = argparse.ArgumentParser(description='Process some arguments.')
parser.add_argument('--file',help='file names to be evaluated',nargs='+')
parser.add_argument('--pot',help='potential at which FED should be shown')
args=parser.parse_args()
args.pot=float(args.pot)

print('%'*20)
print('This script plots the FED corrected for interactions.')
print('%'*20)

def read_output(fname):
    start_reading=False
    energies={}
    for line in open(fname,'r'):
#        if '| CM |  - desc' in line:
        if ' - desc =' in line:
            cpot=line.split('=')[-1].strip()
        if all([a in line for a in ['species','energy']]):
            start_reading=True
            continue
        if start_reading and 'END OUTPUT' in line:
            start_reading=False
            continue
        if start_reading and not '------' in line:
            ls=line.split()
            if cpot not in energies:
                energies[cpot]={}
            energies[cpot][ls[0]]=float(ls[1])
    return energies

def read_gas_output(fname):
    gas_energies={}
    for line in open(fname,'r'):
        ls=line.split()
        if len(ls)<2:
            continue
        if '_g' in ls[0]:
            gas_energies[ls[0]]=float(ls[1])
    return gas_energies

linestyles=cycle(['-','--',':'])

for fname in args.file:
    ls=next(linestyles)
    searchedfile=glob(fname+'/catmap_input/*/*.mkm')
    try:
        mkm_file = sorted( searchedfile, key = lambda file: os.path.getctime(file))[-1]
        model = ReactionModel(setup_file = mkm_file)
    except NameError:
        mkm_file = sorted( searchedfile, key = lambda file: os.path.getctime(file))[-3]
        model = ReactionModel(setup_file = mkm_file)
    
    rxn=model.rxn_expressions
    mec=model.rxn_mechanisms
    
    folder=glob(fname+'/catmap_input/*')
    if fname==args.file[0]:
        min_pot=np.inf
        for f in folder:
            pot=float(f.split('/')[-1].split('_')[-1])
            if abs(pot-args.pot)<min_pot:
                min_pot=abs(pot-args.pot)
                pot_sel=pot
    energies=read_output(fname+'/transport_id000.log')
    gas_energies=read_gas_output(glob(fname+'/catmap_input/*/FED_gas.txt')[0])
    try:
        en_sel=energies[str(pot_sel)]
    except KeyError:
        print("No such key {}. List of keys = {}".format(str(pot_sel),[key for key in energies]))
        sys.exit()
    #ma=MechanismAnalysis(reaction_model=model)
    colors=cycle(['C'+str(i) for i in range(10)])
    nel_count_after=0
    noh_count_after=0
    for m in mec:
        color=next(colors)
        nel_count=0
        noh_count=0
        pos=-1
        pos_ts=0
        en_final_before=0
        for istep,mstep in enumerate(mec[m]):
            nreac=len(rxn[mstep-1].split('->'))
            for ir,r in enumerate(rxn[mstep-1].split('->')):
                itmp=[]
                energy=0
                for rr in r.split('+'):
                    tmp_long=rr.strip().strip('<').split(';')[0]
                    tmp=re.findall('(?:[0-9])?([a-zA-Z-_0-9\*]+)',tmp_long)[0]
                    tmp=tmp.replace('*','')
                    itmp.append(tmp)
                #itmp is now the list of intermediates of this step
                    
                    if tmp in en_sel:
                        energy+=en_sel[tmp]
                    elif tmp in gas_energies: # and tmp not in ['OH_g','ele_g']:
                        #get number count
                        sol=re.findall('([0-9]{0,20})[ ]?'+tmp,tmp_long)
                        if len(sol)>0 and len(sol[0])>0:
                            count=int(sol[0])
                        else:
                            count=1
                        energy+=gas_energies[tmp]*count

                main_reactant=[a for a in itmp if\
                    (not len(re.findall('^[0-9]{0,10}_t$',a))>0) and\
                    (not len(re.findall('^[0-9]{0,10}OH_g$',a))>0) and\
                    (not len(re.findall('^[0-9]{0,10}H_g$',a))>0) and\
                    (not len(re.findall('^[0-9]{0,10}ele_g$',a))>0) and\
                    (not len(re.findall('^[0-9]{0,10}H2O_g$',a))>0)]
                if len(main_reactant)>0:
                    main_reactant='+'.join(main_reactant)
                elif len(main_reactant)==0:
                    main_reactant=''
    
                ##
                #istep counts the reaction step number
                #ir counts the state - IS, TS and FS
                ##
                if istep>0 and ir==0:
                    #save the final energy level from the step before
                    en_final_before=en_final
                if istep==0 and ir==0:
                    en_init=energy #+gas_energies['ele_g']*nel_count_before+gas_energies['OH_g']*noh_count_before
                    pos+=1
                    offset=en_init
                    plt.plot([pos,pos+0.5],[0,0], marker='None', ls=ls, lw=2.5, c=color)
                    plt.annotate(xy=[pos+0.25,0.1], s=main_reactant,color=color,ha='center')
                elif istep>0 and ir==0:
                    en_init=energy
                elif nreac==3 and r==rxn[mstep-1].split('->')[1]:
                    pos_ts=pos+0.5
                    en_ts=(energy-en_init)+en_final_before #+gas_energies['ele_g']*nel_count_before+gas_energies['OH_g']*noh_count_before
                if r==rxn[mstep-1].split('->')[-1]:
                    en_final=(energy-en_init)+en_final_before #+gas_energies['ele_g']*nel_count+gas_energies['OH_g']*noh_count
                    pos+=1
                    plt.plot([pos,pos+0.5],[en_final,en_final], marker='None', ls=ls, lw=2.5, c=color)
                    plt.annotate(xy=[pos+0.25,en_final+0.1], s=main_reactant,color=color,ha='center')
    
            if nreac==3 and en_ts>max(en_final_before,en_final): # and r==rxn[mstep-1].split('->')[1]:
                #transition state
                x=np.linspace(pos_ts,pos_ts+0.5,100)
                a=(en_final_before+en_final-2*en_ts)/0.25**2/2.
                b=(en_ts-en_final_before+a*0.25**2)/0.25
                c=en_ts
                plt.plot(x,a*(x-(pos_ts+0.25))**2+b*(x-(pos_ts+0.25))+c,ls,lw=1.5,c=color)
            else:
                pos_ts=pos-0.5
                plt.plot([pos_ts,pos_ts+0.5],[en_final_before,en_final],ls,lw=1.5,c=color)
plt.title('FED at {} V vs. SHE'.format(pot_sel))
plt.ylim([-8,1])
plt.show()
