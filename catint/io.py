import pickle
from shutil import copyfile as copy

def save_obj(folder,obj, name ):
    with open(folder+'/'+ name + '.pkl', 'wb') as f:
        pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)

def load_obj(name,fname):
    with open(fname+'/' + name + '.pkl', 'rb') as f:
        return pickle.load(f)

def replace_line(file_name, line_num, text):
    lines = open(file_name, 'r').readlines()
    lines[line_num] = text
    out = open(file_name, 'w')
    out.writelines(lines)
    out.close()

def insert_line(file_name,line_num,text):
    f = open(file_name, "r")
    contents = f.readlines()
    f.close()
    contents.insert(line_num, text)
    f = open(file_name, "w")
    contents = "".join(contents)
    f.write(contents)
    f.close() 

def save_all(tp,only=None):
    """saves all dictionaries and arrays into pickle files in the obj folder"""
    if only is None:
        save_obj(tp.outputfoldername,tp.all_data,'all_data')
        save_obj(tp.outputfoldername,tp.species,'species')
        save_obj(tp.outputfoldername,tp.system,'system')
        save_obj(tp.outputfoldername,tp.descriptors,'descriptors')
        save_obj(tp.outputfoldername,tp.xmesh,'xmesh')
        save_obj(tp.outputfoldername,tp.tmesh,'tmesh')
        save_obj(tp.outputfoldername,tp.electrode_reactions,'electrode_reactions')
        save_obj(tp.outputfoldername,tp.electrolyte_reactions,'electrolyte_reactions')
        save_obj(tp.outputfoldername,tp.comsol_outputs,'comsol_outputs')
        save_obj(tp.outputfoldername,tp.comsol_outputs_data,'comsol_outputs_data')
    else:
        save_obj(tp.outputfoldername,tp.all_data,only)

def add_to_dict(self,append,*expr):
    (filename,line_number,function_name,text)=traceback.extract_stack()[-2]
    begin=text.find('add_to_dict(')+len('add_to_dict(')
    end=text.find(')',begin)
    text=[name.strip()+append for name in text[begin:end].split(',')]
    text=text[1:]
    data_dict=dict(zip(text,expr))
    if not hasattr(self,'data_dict'):
        self.data_dict={}
    self.data_dict.update(data_dict)
    for key in self.data_dict:
        if not hasattr(self,key):
            setattr(self,key,self.data_dict[key])
    return self.data_dict

def read_all(tp,fname,only=None):
    if only is None:
        tp.all_data=load_obj('all_data',fname)
        tp.species=load_obj('species',fname)
        tp.system=load_obj('system',fname)
        tp.descriptors=load_obj('descriptors',fname)
        tp.xmesh=load_obj('xmesh',fname)
        tp.tmesh=load_obj('tmesh',fname)
        tp.electrode_reactions=load_obj('electrode_reactions',fname)
        tp.electrolyte_reactions=load_obj('electrolyte_reactions',fname)
        tp.xmax=max(tp.xmesh)
        tp.nx=len(tp.xmesh)
        tp.dx=tp.xmesh[1]-tp.xmesh[0]
        tp.tmax=max(tp.tmesh)
        tp.nt=len(tp.tmesh)
        tp.dt=tp.tmesh[1]-tp.tmesh[0]
    else:
        if type(only)==list:
            for o in only:
                value=load_obj(o,fname)
                setattr(tp,o,value)
        elif type(only)==str:
            value=load_obj(only,fname)
            setattr(tp,only,value)

def print_diff(dictname,varname,missing):
    print('  DIFF -- {} key of {} is missing in model {}'.format(varname,dictname,missing))

#def diff(tp1,tp2):
#    """shows difference between two transport models"""
#    print 'Evaluating differences between transport models.'
#    for a in tp1.all_data:
#        if a not in tp2.all_data:
#            print_diff('all_data',a,'tp2')
#        else:
#            for b in tp1.all_data[a]:
#                
#        for b in tp2.all_data[a]:
#
