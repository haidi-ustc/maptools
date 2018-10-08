#!/usr/bin/env python
from   __future__ import unicode_literals, print_function
import os
import json
import itertools
from pymatgen.symmetry import analyzer
from six.moves import input
from pymatgen.io.vasp.sets import *
from pymatgen.io.vasp.inputs import Poscar, Kpoints,Potcar
from pymatgen import Specie, Lattice, Structure,Molecule
from pymatgen.io.vasp.outputs import Vasprun
from pymatgen.electronic_structure.plotter import DosPlotter,BSPlotter
from maptool.read_structure import readstructure
from maptool.utils import sepline,wait_sep,procs
from maptool.vaspout import *
from maptool.describe_file import *
from maptool.molecule_analysis import MoleculeAnalysis
from maptool.crystal_analysis import CrystalAnalysis
from maptool.maptool_configuration import MMolecule,MStructure
from maptool.mathematics import  unit_vector
from maptool import NAME

def is_pbc():
    print('your choice ?')
    print('{} >>> {}'.format('1','crystal'))
    print('{} >>> {}'.format('2','molecule/cluster'))
    wait_sep()
    in_str=""
    while in_str=="":
       in_str=input().strip()
    choice=int(in_str)
    assert choice in [1,2]
    if choice==1:
       return True
    else:
       return False

def structure_symmetry():
    if is_pbc():
       struct=readstructure(crystal=True,molecule=False)
       ast=analyzer.SpacegroupAnalyzer(struct)
       print("{} : {}".format('Structure Type','periodicity'))
       print("{} : {}".format('Lattice Type',ast.get_lattice_type()))
       print("{} : {}".format('Space Group ID',ast.get_space_group_number()))
       print("{} : {}".format('International Symbol',ast.get_space_group_symbol()))
       print("{} : {}".format('Hall Symbol',ast.get_hall()))
       return
    else:
       struct=readstructure(crystal=False,molecule=True)
       ast=analyzer.PointGroupAnalyzer(struct)
       print("{} : {}".format('Structure Type','non-periodicity'))
       print("{} : {}".format('International Symbol',ast.get_pointgroup()))
       return

def structure_finger_print():
    h_str="# "
    if is_pbc():
       struct=readstructure(crystal=True,molecule=False,cano=True)
       struct_ana=CrystalAnalysis(struct)
       (rdf_x,rdf_species)=struct_ana.fp_oganov()
    else:
       struct=readstructure(crystal=False,molecule=True,cano=True)
       struct_ana=MoleculeAnalysis(struct)    
       (rdf_x,rdf_species)=struct_ana.discrete_radial_distribution_function()

    data=np.zeros((len(rdf_x),len(rdf_species)+1))
    data[:,0]=rdf_x
    tmp1_str="#%(key1)+12s"
    tmp2_dic={'key1':'Distance/A'}
    for i,key in enumerate(rdf_species.keys()):
        data[:,i+1]=rdf_species[key]
        tmp1_str+="%(key"+str(i+2)+")+12s"
        if (isinstance(struct,MStructure)):
           tmp2_dic["key"+str(i+2)]=str(key)
        if (isinstance(struct,Molecule)):
           tmp2_dic["key"+str(i+2)]=str((key[0].number,key[1].number))

    for el in struct.types_of_specie:
        h_str+=str(el)+' |->'+str(el.number)+' '

    head_line=h_str+'\n'+tmp1_str % tmp2_dic

    filename="FingerPrint.dat"
    proc_str="Writting data to "+ filename +" File ..."
    procs(proc_str,0,sp='-->>') 
    write_col_data(filename,data,head_line=head_line)

def structures_difference(distance_tolerance=0.1,rcut=30):
    if is_pbc():
       structs,fnames=readstructure(crystal=True,molecule=False,multi_files=True,cano=True)
    else:
       structs,fnames=readstructure(crystal=False,molecule=True,multi_files=True,cano=True) 
   
    for i in structs:
        print(i)
    for i in fnames:
        print(i)

    n_struct=len(structs)
    proc_str="Total Number of Files are : "+str(n_struct)
    procs(proc_str,0,sp='-->>')
    
    dij=np.zeros((n_struct,n_struct))

    for i, j in itertools.combinations(range(n_struct), 2):
        struct_i=structs[i]
        struct_j=structs[j]
        if isinstance(struct_i,MStructure):
           dij[i,j]=distance(struct_i,struct_j,rcut,pbc=True)
           dij[j,i]=dij[i,j]
        else:
           dij[i,j]=distance(struct_i,struct_j,rcut,pbc=False)
           dij[j,i]=dij[i,j]
    filename="StructureDistance.dat"
    proc_str="Writting data to "+ filename +" File ..."
    procs(proc_str,0,sp='-->>') 
    write_col_data(filename,dij,head_line='')

    data=''
    for i,st_name in enumerate(fnames):
        data+="%(key1)+5s  %(key2)+s"%{'key1':str(i),'key2':fnames[i]}+'\n'
    filename="StructureList.dat"
    proc_str="Writting data to "+ filename +" File ..."
    procs(proc_str,0,sp='-->>') 
    head_line="#%(key1)+5s %(key2)+s"%{'key1':'index','key2':'structure_name'}+'\n'
    data=head_line+data
    write_col_data(filename,data,str_data=True)
       
def distance(struct1,struct2,rcut,pbc=False,):

    fingerprints = {}
    for struct in [struct1,struct2]:
       if pbc:
          analysis = CrystalAnalysis(struct, radius=rcut)
       else:
          analysis = MoleculeAnalysis(struct)
           
       x, ys = analysis.fp_oganov()
       fingerprint = {'id': struct}
       for k in ys:
           atomic_number1 = k[0]
           atomic_number2 = k[1]
           pair = '%06d' % min(atomic_number1 * 1000 + atomic_number2,atomic_number2 * 1000 + atomic_number1)
           fingerprint[pair] = list(ys[k])
       fingerprints[id(struct)] = fingerprint

    dij = []
    for pair in fingerprints[id(struct1)]:
        if pair in fingerprints[id(struct2)] and pair != 'id':
            uvect1 = unit_vector(fingerprints[id(struct1)][pair])
            len1=len(uvect1)
            uvect2 = unit_vector(fingerprints[id(struct2)][pair])
            len2=len(uvect2)
            if len1>len2:
               vect=np.zeros(len1) 
               vect[0:len2]=uvect2
               uvect2=vect
            elif len1<len2:
               vect=np.zeros(len2) 
               vect[0:len1]=uvect1
               uvect1=vect
            else:
               pass
            dij.append(0.5 * (1.0 - np.dot(uvect1, uvect2)))
    distance = float(np.mean(dij))
    return distance

def get_primitive_cell():
    struct=readstructure()
    sepline(ch='Primitive Cell',sp='-')
    ast=analyzer.SpacegroupAnalyzer(struct)
    prim_st=ast.find_primitive()
    print(prim_st)
    sepline()
    print('save to '+NAME+'_primitive.vasp')
    prim_st.to(filename=NAME+'_primitive.vasp',fmt='poscar')
    sepline()
    return

def get_conventional_cell():
    struct=readstructure()
    sepline(ch='conventional  cell',sp='-')
    ast=analyzer.SpacegroupAnalyzer(struct)
    conv_st=ast.get_conventional_standard_structure()
    print(conv_st)
    sepline()
    print('save to '+NAME+'_convention.vasp')
    conv_st.to(filename=NAME+'_conventional.vasp',fmt='poscar')
    sepline()
    return 
