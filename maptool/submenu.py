#!/usr/bin/env python
from __future__ import unicode_literals, print_function
from six.moves import input
from maptool.operation import    random_operation, covert_operation,cleave_operation, strain_operation,\
                                 twoD_operation, build_operation,get_xrd
from maptool.online_extract import  online_get_banddos, online_get_structure, online_get_properties
from maptool.read_structure import readstructure
from maptool.utils import sepline,wait_sep,procs
from maptool.describe_file import *
from maptool.vaspinput import generate_all_input,generate_incar,generate_kpoint,generate_potcar
from maptool.vaspout   import total_dos, projected_dos, band_structure, projected_band_structure, select_one_band_structure,\
                              charge_density, spin_density, charge_density_diff, spin_density_component, chg_locp_average,\
                              optics_analysis, aimd_analysis, elastic_analysis
from maptool.function  import structure_symmetry,structure_finger_print,structures_difference,distance,get_primitive_cell,get_conventional_cell

def select_function(choice):

# structure operation
    if choice=="a1":
       random_operation()
    elif choice=="a2":
       covert_operation()
    elif choice=="a3":
       build_operation()
    elif choice=="a4":
       cleave_operation()
    elif choice=="a5":
       strain_operation()
    elif choice=='a6':
       twoD_operation()

# structure analysis
    elif choice=="b1":
       structure_symmetry()
    elif choice=="b2":
       structure_finger_print()
    elif choice=="b3":
       structures_difference()   
    elif choice=="b4":
       get_primitive_cell()
    elif choice=="b5":
       get_conventional_cell()
    elif choice=="b6":
       get_xrd()

# vasp in/out tools
    elif choice=="c1":
       struct=readstructure(crystal=True,molecule=False,filename='POSCAR',cano=False)
       sepline(ch=' prepare intput files ',sp='-')
       print('your choce ?')
       print('{} >>> {}'.format('1','prepare all files automatically'))
       print('{} >>> {}'.format('2','prepare INCAR file'))
       print('{} >>> {}'.format('3','prepare KPOINTS file'))
       print('{} >>> {}'.format('4','prepare POTCAR file'))
       wait_sep()
       in_str=""
       while in_str=="":
          in_str=input().strip()
       choice=int(in_str)
       if choice==1:
          generate_all_input(struct)
       elif choice==2:
          generate_incar(struct)
       elif choice==3:
          generate_kpoint(struct)
       elif choice==4:
          generate_potcar(struct)
       else:
          print("unknown choice")
    
    elif choice=="cxx":
       sepline(ch=' summary output files ',sp='=')
       print('{} >>> {}'.format('1','describe OUCAR file'))
       print('{} >>> {}'.format('2','describe OSICAR file'))
       print('{} >>> {}'.format('3','describe vasprun.xml file'))
       wait_sep()
       in_str=""
       while in_str=="":
          in_str=input().strip()
       choice=int(in_str)
       if choice==1:
          describe_outcar()
       elif choice==2:
          describe_OSICAR()
       elif choice==3:
          describe_vasprun()
       else:
          print("unknow choice")

    elif choice=="c2":
       sepline(ch=' vasp output analysis ',sp='-')
       print('{} >>> {}'.format('1 ','total density of states'))
       print('{} >>> {}'.format('2 ','projected density of states'))
       print('{} >>> {}'.format('3 ','band structure'))
       print('{} >>> {}'.format('4 ','projected band structure'))
       print('{} >>> {}'.format('5 ','select one band structure'))
       print('{} >>> {}'.format('6 ','charge density'))
       print('{} >>> {}'.format('7 ','spin density'))
       print('{} >>> {}'.format('8 ','charge density difference'))
       print('{} >>> {}'.format('9 ','spin density component: up/down'))
       print('{} >>> {}'.format('10','average charge density/potential'))
       print('{} >>> {}'.format('11','optics analysis'))
       print('{} >>> {}'.format('12','mechanical analysis'))
       print('{} >>> {}'.format('13','ab initio molecular dynamics analysis'))
       wait_sep()
       in_str=""
       while in_str=="":
          in_str=input().strip()
       choice=int(in_str)
       if choice==1:
          total_dos()
       elif choice==2:
          projected_dos()
       elif choice==3:
          band_structure()
       elif choice==4:
          projected_band_structure()
       elif choice==5:
          select_one_band_structure()
       elif choice==6:
          charge_density()
       elif choice==7:
          spin_density()
       elif choice==8:
          charge_density_diff()
       elif choice==9:
          spin_density_component()
       elif choice==10:
          chg_locp_average()
       elif choice==11:
          optics_analysis()
       elif choice==12:
          elastic_analysis()
       elif choice==13:
          aimd_analysis()
       else:
          print("unknow choice")

# online exctraction
    elif choice=="e1":
       online_get_banddos()
    elif choice=="e2":
       online_get_structure()
    elif choice=="e3":
       online_get_properties()

    elif choice=="88":
         return
    else:
       print("unknow choice")
       return

