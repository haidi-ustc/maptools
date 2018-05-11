#!/usr/bin/env python
from six.moves import input
import uuid
import os
import glob
from maptool import mpt_log,NAME
from maptool.utils import sepline,wait_sep
from pymatgen import Structure,Molecule
from maptool.maptool_configuration import MStructure,MMolecule

infile_list_com=[NAME+'.nc',NAME+'.yaml',NAME+'.json']
infile_list_c=['POSCAR','CONTCAR','CHGCAR','LOCPOT',NAME+'.cif',NAME+'.xsf']+infile_list_com
infile_list_m=[NAME+'.xyz',NAME+'.mol' ]+infile_list_com

def readstructure(crystal=True,molecule=False,filename=None,multi_files=False,cano=False):

    if multi_files:
       if crystal:
          print('input the full file names seperated by space (only POSCAR format)')
          print('supported format 1: a.vasp b.vasp')
          print('supported format 2: *.vasp')
          structs=[] 
          wait_sep()
          in_str=""
          while in_str=="":
                in_str=input().strip()
          if '*' in in_str:
              fnames=glob.glob(in_str)
          else:
              fnames=in_str.split()

          for fname in fnames:
              nfname=str(uuid.uuid4())+'_POSCAR'
              os.symlink(fname,nfname)
              mpt_log.debug("linke file for %s is %s" % (fname,nfname))
              fname=nfname
              struct=MStructure.from_file(fname)
              if cano:
                 struct.canonical_form()
              os.unlink(fname)
              structs.append(struct)
          return structs,fnames
          
       if molecule:
          print('input the full file names seperated by space (only for .xyz format)')
          print('supported format 1: a.xyz b.xyz')
          print('supported format 2: *.xyz')
          structs=[]
          wait_sep()
          in_str=""
          while in_str=="":
                in_str=input().strip()
          if '*' in in_str:
              fnames=glob.glob(in_str)
          else:
              fnames=in_str.split()

          for fname in fnames:
              struct=MMolecule.from_file(fname)
              if cano:
                 struct.canonical_form()
              structs.append(struct)
          return structs,fnames

    if filename is not None:
       try:
          struct=MStructure.from_file(filename)
          if cano:
              struct.canonical_form()
       except:
          struct=MMolecule.from_file(filename)
          if cano:
              struct.canonical_form()
       return struct

#    sepline(ch='Read Structure',sp='=')
    print('your choice ?')
    print('{} >>> {}'.format('0','specific the file name!'))
    if crystal:
       print('{} >>> {}'.format('1','POSCAR'))
       print('{} >>> {}'.format('2','CONTCAR'))
       print('{} >>> {}'.format('3','CHGCAR'))
       print('{} >>> {}'.format('4','LOCPOT'))
       print('{} >>> {}'.format('5',NAME+'.cif'))
       print('{} >>> {}'.format('6',NAME+'.xsf'))
       print('{} >>> {}'.format('7',NAME+'.nc'))
       print('{} >>> {}'.format('8',NAME+'.yaml'))
       print('{} >>> {}'.format('9',NAME+'.json'))
    if molecule:
       print('{} >>> {}'.format('1',NAME+'.xyz'))
       print('{} >>> {}'.format('2',NAME+'.mol'))
       print('{} >>> {}'.format('3',NAME+'.nc'))
       print('{} >>> {}'.format('4',NAME+'.yaml'))
       print('{} >>> {}'.format('5',NAME+'.json'))

    wait_sep()
    in_str=""
    while in_str=="":
          in_str=input().strip()
    fstin=int(in_str)
    if crystal:
       if fstin==0:
          print('input the full file name')
          wait_sep()
          in_str=""
          while in_str=="":
                in_str=input().strip()
          fname=in_str
          if fname.endswith('.vasp'):
             nfname=str(uuid.uuid4())+'_POSCAR'
             os.symlink(fname,nfname)
             fname=nfname 
       else:
          fname=infile_list_c[fstin-1]
       struct=MStructure.from_file(fname)
       if cano:
          struct.canonical_form()

    if molecule:
        if fstin==0:
           print('input the full file name')
           wait_sep()
           in_str=""
           while in_str=="":
                in_str=input().strip()
           fname=in_str
        else:
           fname=infile_list_m[fstin-1]
        struct=MMolecule.from_file(fname)
        if cano:
           struct.canonical_form()
    return struct

if __name__=='__main__':
   structs=readstructure(crystal=True,molecule=False,filename=None,multi_files=True)
