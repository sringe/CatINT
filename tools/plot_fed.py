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

searchedfile=glob(args.file[0]+'/catmap_input/*/*.mkm')
try:
    mkm_file = sorted( searchedfile, key = lambda file: os.path.getctime(file))[-1]
    model = ReactionModel(setup_file = mkm_file)
except NameError:
    mkm_file = sorted( searchedfile, key = lambda file: os.path.getctime(file))[-3]
    model = ReactionModel(setup_file = mkm_file)

rxn=model.rxn_expressions
mec=model.rxn_mechanisms

folder=glob(args.file[0]+'/catmap_input/*')
min_pot=np.inf
for f in folder:
    pot=float(f.split('/')[-1].split('_')[-1])
    if abs(pot-args.pot)<min_pot:
        min_pot=abs(pot-args.pot)
        pot_sel=pot
#pH=model.bulk_ph
pH=7
energies=read_output(args.file[0]+'/transport_id000.log')
try:
    en_sel=energies[str(pot_sel)]
except KeyError:
    print("No such key {}. List of keys = {}".format(str(pot_sel),[key for key in energies]))
    sys.exit()
#ma=MechanismAnalysis(reaction_model=model)
colors=cycle(['C'+str(i) for i in range(10)])
for m in mec:
    color=next(colors)
    print m
    nel_count=0
    noh_count=0
    pos=-1
    pos_ts=0
    for istep,mstep in enumerate(mec[m]):
        nreac=len(rxn[mstep-1].split('->'))
        for ir,r in enumerate(rxn[mstep-1].split('->')):
            itmp=[]
            energy=0
            for rr in r.split('+'):
                tmp=''.join(re.findall('([0-9a-zA-Z-_]+)',rr.strip().strip('<').split(';')[0]))
                itmp.append(tmp)
            #itmp is now the list of intermediates of this step
                if tmp in en_sel:
                    energy+=en_sel[tmp]
            nel=0
            noh=0
            for tmp in itmp:
                if 'ele_g' in tmp:
                    sol=re.findall('([0-9]+)',tmp)
                    if len(sol)>0:
                        nel=int(sol[0])
                    else:
                        nel=1
                if 'OH_g' in tmp:
                    sol=re.findall('([0-9]+)',tmp)
                    if len(sol)>0:
                        noh=int(sol[0])
                    else:
                        noh=1
            nel_count_before=nel_count
            nel_count+=nel
            noh_count_before=noh_count
            noh_count+=noh
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

            if istep>0 and ir==0:
                en_init=en_final
            if istep==0 and ir==0:
                en_init=energy+pot_sel*nel_count_before+noh_count_before*0.059*pH
                pos+=1
                plt.plot([pos,pos+0.5],[en_init,en_init], marker='None', ls='-', lw=2.5, c=color)
                plt.annotate(xy=[pos+0.25,en_init+0.1], s=main_reactant,color=color,ha='center')
            elif nreac==3 and r==rxn[mstep-1].split('->')[1]:
                pos_ts=pos+0.5
                en_ts=energy+pot_sel*nel_count_before+noh_count_before*0.059*pH
            if r==rxn[mstep-1].split('->')[-1]:
                en_final=energy+pot_sel*nel_count+noh_count*0.059*pH
                pos+=1
                plt.plot([pos,pos+0.5],[en_final,en_final], marker='None', ls='-', lw=2.5, c=color)
                plt.annotate(xy=[pos+0.25,en_final+0.1], s=main_reactant,color=color,ha='center')

        print r,en_init,en_ts,en_final,'pot=',pot_sel
        if nreac==3 and en_ts>max(en_init,en_final): # and r==rxn[mstep-1].split('->')[1]:
            #transition state
            x=np.linspace(pos_ts,pos_ts+0.5,100)
            a=(en_init+en_final-2*en_ts)/0.25**2/2.
            b=(en_ts-en_init+a*0.25**2)/0.25
            c=en_ts
            plt.plot(x,a*(x-(pos_ts+0.25))**2+b*(x-(pos_ts+0.25))+c,'-',lw=1.5,c=color)
        else:
            pos_ts=pos-0.5
            plt.plot([pos_ts,pos_ts+0.5],[en_init,en_final],'-',lw=1.5,c=color)
plt.ylim([-8,1])
plt.show()
