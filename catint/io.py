import pickle

def save_obj(obj, name ):
    with open('obj/'+ name + '.pkl', 'wb') as f:
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

def save_all(tp):
    """saves all dictionaries and arrays into pickle files in the obj folder"""
    save_obj(tp.all_data,'all_data')
    save_obj(tp.species,'species')
    save_obj(tp.system,'system')
    save_obj(tp.descriptors,'descriptors')
    save_obj(tp.xmesh,'xmesh')
    save_obj(tp.tmesh,'tmesh')
    save_obj(tp.electrode_reactions,'electrode_reactions')
    save_obj(tp.electrolyte_reactions,'electrolyte_reactions')

def read_all(tp,fname):
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

