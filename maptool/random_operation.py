#!/usr/bin/env python
from six.moves import input
import random
import numpy as np
from pymatgen import Structure,Molecule,Specie,Composition,Lattice
from maptool.read_structure import readstructure
from maptool.utils import *
from maptool.mathematics import apply_rotation
from maptool.function import is_pbc

class RandomStructure(Structure):

    def __init__(self,composition,method='scaling',stabiliztion_number=3,paral=5,periodic=True,factor_opptimal_volume=8):
    
        self._comp=Composition(composition)
        self._method=method
        self.method_list=['scaling','stretching','pseudorandom']
    
    @property
    def composition(self):
        return self._comp
    
    @property
    def method(self):
        return self._method

    @staticmethod
    def get_one_random_structure():
        pass

    def get_random_structures(number):
        structs=[]
        for i in range(number):
            structs.append(get_one_random_structure())
        return structs

class RandomMolecule(Molecule):

    def __init__(self,composition,method='scaling',stabiliztion_number=3,paral=5,periodic=True,factor_opptimal_volume=8):
    
        self._comp=Composition(composition)
        self._method=method
        self.method_list=['scaling','stretching','pseudorandom']
    
    @property
    def composition(self):
        return self._comp
    
    @property
    def method(self):
        return self._method

    @staticmethod
    def get_one_random_molecule():
        pass

    def get_random_molecules(number):
        molecules=[]
        for i in range(number):
            molecules.append(get_one_random_structure())
        return molecules

def random_purturbation_index():
    if is_pbc():
       struct=readstructure(crystal=True,molecule=False)
       natom=len(struct)
       d=np.zeros(2,dtype=int)
       while d[0]==d[1]:
           d=np.random.randint(0,natom,2)
       coord=struct.frac_coords
       coord[[d[0], d[1]], :] = coord[[d[1], d[0]], :]
       tmp_struct=Structure(struct.lattice,struct.species,coord)
       fname='swap_'+str(d[0]+1)+'_'+str(d[1]+1)+'.vasp'
       proc_str="Saving data to "+ fname +" File ..."
       procs(proc_str,0,sp='-->>')
       tmp_struct.to(filename=fname,fmt='poscar')
    else:
       struct=readstructure(crystal=False,molecule=True)
      # print(struct)     
       coord=struct.cart_coords
       natom=len(struct)
       d=np.zeros(2,dtype=int)
       while d[0]==d[1]:
           d=np.random.randint(0,natom,2)
       coord[[d[0], d[1]], :] = coord[[d[1], d[0]], :]
       tmp_struct=Molecule(struct.species,coord)
       fname='swap_'+str(d[0]+1)+'_'+str(d[1]+1)+'.xyz'
       proc_str="Saving data to "+ fname +" File ..."
       procs(proc_str,0,sp='-->>')
       tmp_struct.to(filename=fname,fmt='xyz')
    return

def random_disturbing_lat():
    struct=readstructure(crystal=True,molecule=False)
    print("input the maximum displacement and with fix diagonal ")
    print("or non-diagonal element like this:  0.01 T F")
    print("it means maximum displacement is 0.01 ")
    print("the diagonal element will be fixed,") 
    print("while random disturbing will be add to non-digaonal element")
    wait_sep()
    in_str=""
    while in_str=="":
       in_str=input().strip().split()
#    print(in_str)
    maxdelta=float(in_str[0])
    diag=in_str[1].lower()=='t'
    nondiag=in_str[2].lower()=='t'
#    print(diag)
#    print(nondiag)
    stress_eps = np.random.random(6) * 2 * maxdelta - maxdelta
    if diag:
       stress_eps[:3] = 0
    if nondiag:
       stress_eps[-3:] = 0
    new_lattice=deform_cell(struct,stress_eps)
    tmp_struct=Structure(new_lattice,struct.species,struct.frac_coords)
    fname='random.vasp'
    proc_str="Saving data to "+ fname +" File ..."
    procs(proc_str,0,sp='-->>')
    tmp_struct.to(filename=fname,fmt='poscar')
    return

def deform_cell(struct, stress_eps):
    stress = np.eye(3) + np.diag(stress_eps[:3]) + np.array([[0.0, stress_eps[3], stress_eps[4]],
                                                             [stress_eps[3], 0.0, stress_eps[5]],
                                                             [stress_eps[4], stress_eps[5], 0.0]])
    return np.dot(stress, struct.lattice.matrix)

def random_disturbing_pos():

    def random_move_one_atom(coords, mu=0.1, sigma=0.01):

        index = random.randint(0, len(coords) - 1)
        radius = np.abs(np.random.normal(mu, sigma))
        theta_x = 2 * np.pi * np.random.random_sample()
        theta_y = 2 * np.pi * np.random.random_sample()
        theta_z = 2 * np.pi * np.random.random_sample()
        vector = apply_rotation([1, 0, 0], theta_x, theta_y, theta_z)
        coords[index] += vector*radius
        return coords
    
    print("input the maximum displacement(<0.25 in Angstrom)")
    wait_sep()
    in_str=""
    while in_str=="":
       in_str=input().strip()
    epsilon=float(in_str)
    assert epsilon < 0.3
    
    if is_pbc():
       struct=readstructure(crystal=True,molecule=False)
       coords=struct.cart_coords
       for iatom in range(len(struct)):
           coords=random_move_one_atom(coords,mu=epsilon)
       tmp_struct=Structure(struct.lattice,struct.species,coords,coords_are_cartesian=True)
       fname='random.vasp'
       proc_str="Saving data to "+ fname +" File ..."
       procs(proc_str,0,sp='-->>')
       tmp_struct.to(filename=fname,fmt='poscar')
    else:
       struct=readstructure(crystal=False,molecule=True)
       coords=struct.cart_coords
       for iatom in range(len(struct)):
           coords=random_move_one_atom(coords,mu=epsilon)

       tmp_struct=Molecule(struct.species,coords)
       fname='random.xyz'
       proc_str="Saving data to "+ fname +" File ..."
       procs(proc_str,0,sp='-->>')
       tmp_struct.to(filename=fname,fmt='xyz')        
    return
