from subprocess import call
import os
from shutil import copyfile,copytree
import itertools
import sys

desc={
    'H2O-CO-ele*_t':{\
            'beta':[0.3,0.5,0.7],\
            'energy':[-0.4,0.0,0.4]
            },
    'CHO-H2O-ele*_t':{\
            'beta':[0.3,0.5,0.7],\
            'energy':[-0.4,0.0,0.4]
            },
    'CH-OH-ele*_t':{\
            'beta':[0.3,0.5,0.7],\
            'energy':[-0.4,0.0,0.4]
            },
    }


#create results folder:
if not os.path.exists('collected_results'):
    os.makedirs('collected_results')


newdict={}
for d in desc:
    for c in desc[d]:
        newdict[d+';'+c]=desc[d][c]
keys,values=zip(*newdict.items())
to_iter=[dict(zip(keys, v)) for v in itertools.product(*values)]

copyfile('catmap_CO2R_template.mkm.bu','tmp_template.mkm')
copyfile('catmap_CO2R_energies.txt.bu','tmp_energies.txt')

os.chdir('collected_results')
root=os.getcwd()

with open(root+'/log.txt','w') as of:
    for i,d in enumerate(to_iter):
        if i>0:
            sys.exit()
        os.chdir(root)
        if not os.path.isdir('results_'+str(i)):
            os.makedirs('results_'+str(i))
        os.chdir('results_'+str(i))
        print os.getcwd()
        copytree(root+'/../data','data')
        copyfile(root+'/../test_tp_new_CO2.py','test_tp_new_CO2.py')
        copyfile(root+'/../read_data.py','read_data.py')
        copyfile(root+'/../tmp_template.mkm','tmp_template.mkm')
        copyfile(root+'/../tmp_energies.txt','tmp_energies.txt')
#        copyfile(root+'/../slac.small','slac.small')
        for r in to_iter[i]:
            mol=r.split(';')[0]
            prop=r.split(';')[1]
            val=to_iter[i][r]
            of.write('modifying {} {} {}\n'.format(mol,prop,val))
            of.flush()
            if prop in ['beta']:
                with open('catmap_CO2R_template.mkm','w') as of_temp:
                    for line in open('tmp_template.mkm','r'):
                        if '-> '+mol+' <->' in line:
                            ls=line.split(';')[0]+'; beta='+str(val)+'\',\n'
                            of_temp.write(ls) #'\''+'CO*_t + H2O_g + ele_g <-> H2O-CO-ele*_t <-> CHO*_t + OH_g'+'; beta='+desc[d][c])
                        else:
                            of_temp.write(line)
                copyfile('catmap_CO2R_template.mkm','tmp_template.mkm')
            elif prop in ['energy']:
                with open('catmap_CO2R_energies.txt','w') as of_en:
                    for line in open('tmp_energies.txt','r'):
                        if mol.replace('*','').replace('_t','') in line:
                            ls=line.split('\t')
                            ls[3]=str(float(ls[3])+val)
                            of_en.write('\t'.join(ls))
                        else:
                            of_en.write(line)
                copyfile('catmap_CO2R_energies.txt','tmp_energies.txt')
        of.write('ready to start\n')
        of.flush()
        call("esp-ver-bsub 19 test_tp_new_CO2.py",shell=True)
#        call(" " + " slac.small", shell=True)
#        call("python" + " test_tp_new_CO2.py", shell=True)
        #os.rename('CO2R_results','results_'+str(i))
