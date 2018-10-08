#!/usr/bin/env python
import os
from six.moves import input
import numpy as np
from pymatgen import  Element, Structure
from maptool.constants import Len

def sepline(ch='-',sp='-'):
    r'''
    seperate the output by '-'
    '''
    print(ch.center(Len,sp))

def box_center(ch='',fill=' ',sp="|"):
    r'''
    put the string at the center of |  |
    '''
    strs=ch.center(Len,fill)
    print(sp+strs[1:len(strs)-1:]+sp)

def warn_tip(idex=0,instr=''):
    if idex==0:
       stc=" Warning "
    else:
       stc=" Tips "
    box_center(ch=stc,fill='-',sp="+")
    for i_str in instr.split("\n"):
        box_center(ch=' '+i_str+' ',fill=' ',sp=" ")
    box_center(ch='-',fill='-',sp="+")

def wait_sep():
    r'''
    waiting for user input
    '''
    print('--------------->>')

def procs(ch,ind,sp='-->>'):
    r'''
    under processing 
    '''
    if ind !=0:
       nstr=sp+ " ("+str(ind)+') '+ch
    else:
       nstr=sp+"     "+ch
    print(nstr)

def check_file(filename):
    r'''
    checking whether exist the specific file
    '''
    if os.path.exists(filename):
       pass
    else:
       print(filename+" file is not found")
       os._exit()

def check_matplotlib():
    r'''
    loading the plot tool: matplot lib
    '''
 
    try:
       import matplotlib
       #matplotlib.use('Agg')
    except:
       print("you have to install matplotlib")
       os._exit()

def atom_selection(struct):
    r'''
    select atoms by three different schemes:
    1. by atomic index
    2. by element symbol
    3. by fractional coordinates range
    '''
    print("")
    print("    input data according to tips")
    tip="""
     select atoms by following ways:
     1. atomic index in POSCAR
        i.e. :  1 2 4-8 10 12-30
        i.e. :  1 2 4 8 10 
     2. atomic label
        i.e. :  Si  O
     3. atomic position
        i.e. :  0 0.5 | 0.2 0.4 | 0.3 0.7
        this means atoms with 0<x<0.5, 
        0.2<y<0.4 and 0.3<z<0.7 will be seleted
        or just specific the z coordinates,
        i.e. :  ||0.3 0.7
        """
    print(tip)
    wait_sep()
    in_str=""
    while in_str=="":
          in_str=input().strip()
 
    def parse_index(in_str):
        atom_index=[]
        tmp_str=in_str.split()
        #tmp_str=in_str.split()[1:]
        for i in tmp_str:
            if '-' not in i:
               atom_index.append(int(i))
            else:
               atom_index.extend(range(int(i.split('-')[0]),int(i.split('-')[1])+1))
        return [i-1 for i in atom_index]
 
    def parse_label(in_str):
        atom_label= [Element(elem) for elem in in_str.split()]
        #atom_label= [Element(elem) for elem in in_str.split()[1:]]
        atom_index=[]
        for i, site in enumerate(struct.sites):
            if site.specie in atom_label:
               atom_index.append(i)
        return atom_index
 
    def parse_range(in_str):
 
        def check_frac(coord,lim):
            con=[False]*3
            for i in range(3):
               con[i]=coord[i]>=lim[i][0] and coord[i]<=lim[i][1]
            if np.all(con):
               return True
            else:
               return False
 
        coord_range={}
        tmp_str=in_str.split();tmp1_str=' '.join(tmp_str);tmp2_str=tmp1_str.split("|")
        #tmp_str=in_str.split()[1:];tmp1_str=' '.join(tmp_str);tmp2_str=tmp1_str.split("|")
        icount=0
        for i in tmp2_str:
 
            if i=='':
              coord_range[icount]=''
            else:
              coord_range[icount]=[float(x) for x in i.split()]
            icount+=1
        for key in coord_range.keys():
            if coord_range[key]=='':
               coord_range[key]=[0,1]
        atom_index=[]
        for i, site in enumerate(struct.sites):
            if check_frac(site.frac_coords,coord_range):
               atom_index.append(i)
        return atom_index
 
    if "|" in in_str.strip():
       atom_index_list=parse_range(in_str)
    else:
       for str in in_str.strip():
           if str.isalpha():
              atom_index_list=parse_label(in_str)
              return atom_index_list,in_str
       atom_index_list=parse_index(in_str)
    #if in_str.strip().startswith('a'):
    #   atom_index_list=parse_index(in_str)
    #elif in_str.strip().startswith('b'):
    #   atom_index_list=parse_label(in_str)
    #elif in_str.strip().startswith('c'):
    #   atom_index_list=parse_range(in_str)
       return atom_index_list,in_str

