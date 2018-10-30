#!/usr/bin/env python
from    __future__ import division, unicode_literals, print_function
import math
import os
import numpy as np
from six.moves import input
from monty.json import jsanitize
from pymatgen.core import Structure,Molecule,Lattice
from pymatgen.analysis.adsorption import *
from pymatgen.core.surface import generate_all_slabs, SlabGenerator
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.io.vasp import Poscar
from pymatgen.analysis.diffraction.xrd import XRDCalculator
from maptool.utils import *
from maptool.substrate_matcher import make_connect,get_subs
from maptool.read_structure import readstructure
from maptool.data_io import write_col_data,json_store
from maptool.random_operation import get_random_structure,random_purturbation_index,\
                                     random_disturbing_pos,random_disturbing_lat
from   maptool.mathematics import smear

def random_operation():
    print('your choice ?')
    print('{} >>> {}'.format('1','random structure generating'))
    print('{} >>> {}'.format('2','random perturbation for atom index'))
    print('{} >>> {}'.format('3','random disturbing for lattice matrix'))
    print('{} >>> {}'.format('4','random disturbing for atom position'))
    wait_sep()
    in_str=""
    while in_str=="":
       in_str=input().strip()
    choice=int(in_str)
    assert choice in [1,2,3,4]

    if choice==1:
       get_random_structure()
    elif choice==2:
       random_purturbation_index()
    elif choice==3:
       random_disturbing_lat()
    else:
       random_disturbing_pos()

def covert_operation():
    print('your choice ?')
    print('{} >>> {}'.format('1','crystal'))
    print('{} >>> {}'.format('2','single molecule'))
    wait_sep()
    in_str=""
    while in_str=="":
       in_str=input().strip()
    choice=int(in_str)
    assert choice in [1,2]
    if choice==1:
       struct=readstructure(crystal=True,molecule=False)
       print('input the target file name\nremeber including the format sufix')
       print("supported format: .vasp .xsf .cif .nc .json .yaml ")
    else:
       struct=readstructure(crystal=False,molecule=True)
       print('input the target file name\nremeber including the format sufix')
       print("supported format: .xyz .mol .nc .json .yaml ")
         
    wait_sep()
    in_str=""
    while in_str=="":
       in_str=input().strip()
    outfile=in_str
    fmt=outfile.split('.')[-1]
    struct.to(filename=outfile,fmt=fmt)

def build_operation():
    print('your choice ?')
    print('{} >>> {}'.format('1','build supercell'))
    print('{} >>> {}'.format('2','build nanotube'))
    print('{} >>> {}'.format('3','build absorption configuration'))
    wait_sep()
    in_str=""
    while in_str=="":
       in_str=input().strip()
    choice=int(in_str)

    assert choice in [1,2,3]
    if choice==1:
        struct=readstructure()
        wait_sep()
        tip=        """
Several options are available:

a. A full 3x3 scaling matrix defining the linear combination
   the old lattice vectors. E.g., 2 1 0  0 1 0  0 0 3
   generates a new structure with lattice vectors a' =
   2a + b, b' = 3b, c' = c where a, b, and c are the lattice
   vectors of the original structure.
b. An sequence of three scaling factors. E.g., 2 1 1
   specifies that the supercell should have dimensions 2a x b x
   c.
c. A number, which simply scales all lattice vectors by the
   same factor.
        """
        print(tip)
        wait_sep()
        in_str=""
        while in_str=="":
           in_str=input().strip()
        choice=in_str
        scaling_list=[int(x) for x in choice.split()]
        print(scaling_list)
        if len(scaling_list)==1:
            scales=scaling_list[0]
        elif len(scaling_list)==3:
            scales=scaling_list
        elif len(scaling_list)==9:
            scales=[scaling_list[0:3],scaling_list[3:6],scaling_list[6:9]]
        print(scales)
        struct_cp=struct.copy()
        struct_cp.make_supercell(scales)
        print(struct_cp)
        struct_cp.to(filename='supercell.vasp',fmt='poscar')
#        if  struct_cp==struct and :
#            print('unreasonable input')
#            build_operation()
#        else:
#            return
    elif choice==2:
        struct=readstructure()
        wait_sep()
        tip=        """
unsupported now !
                    """
        print(tip)
        return 
    else:
        data={'max_index': 2, 'min_vacum': 20, 'min_slab': 8, 'repeat': [3, 3, 1]}
        def read_adsorb_config(filename):
            with open(filename,'r') as f:
                 datas=f.readlines()
            list_data=[]
            for i in range(len(datas)):
                list_data.append(datas[i][0:datas[i].find('#')].strip().split('='))

            defined_keys=['method','crystal','molecule','max_index','min_vacum','min_slab','repeat']
            data_dict={}
            for key in defined_keys:
                for li in list_data:
                    if key in li[0]:
                       data_dict[key]=li[1]

            data_dict['method']=int(data_dict.get('method').strip())
            data_dict['crystal']=data_dict.get('crystal').strip()
            data_dict['molecule']=data_dict.get('molecule').strip()
            data_dict['max_index']=int(data_dict.get('max_index','1').strip())
            data_dict['min_vacum']=int(data_dict.get('min_vacum','15').strip())
            data_dict['min_slab']=int(data_dict.get('min_slab','5').strip())
            data_dict['repeat']=[int(x) for x in data_dict.get('repeat','1 1 1').strip().split()]
            return data_dict

        def proc_adsorb(cryst,mol,data):
           if data['method'] ==1:
              asf_slab=AdsorbateSiteFinder(cryst)
              ads_sites=asf_slab.find_adsorption_sites()
              ads_structs=asf_slab.generate_adsorption_structures(mol,repeat=data['repeat'])
              for i in range(len(ads_structs)):
                  ads_struct=ads_structs[i]
                  try:
                     miller_str=[str(j) for j in cryst.miller_index]
                  except:
                     miller_str=['adsorb']
                  filename='_'.join(miller_str)+'-'+str(i)+'.vasp'
                  ads_struct.to(filename=filename,fmt='POSCAR')
           else:
              slabs=generate_all_slabs(cryst,max_index=data['max_index'],min_slab_size=data['min_slab'],
                                   min_vacuum_size=data['min_vacum'],lll_reduce=True)
              for slab in slabs:
                  asf_slab=AdsorbateSiteFinder(slab)
                  ads_sites=asf_slab.find_adsorption_sites()
                  ads_structs=asf_slab.generate_adsorption_structures(mol,repeat=data['repeat'])
                  for i in range(len(ads_structs)):
                      ads_struct=ads_structs[i]
                      miller_str=[str(j) for j in slab.miller_index]
                      filename='adsorb'+'_'.join(miller_str)+'-'+str(i)+'.vasp'
                      ads_struct.to(filename=filename,fmt='POSCAR')
 
        filename='adsorb.cfg' 
        if os.path.exists(filename):
           data=read_adsorb_config(filename)
           assert data['method'] in [1,2]
           cryst=readstructure(filename=data['crystal'])
           mol=readstructure(filename=data['molecule'])
           proc_adsorb(cryst,mol,data)
        else:
           print('your choice ?')
           print('{} >>> {}'.format('1','read slab from file'))
           print('{} >>> {}'.format('2','build slab by bulk'))
           wait_sep()
           in_str=""
           while in_str=="":
              in_str=input().strip()
           choice=int(in_str)
           assert choice in [1,2]
           data['method']=choice
           cryst=readstructure(crystal=True,molecule=False)
           mol=readstructure(crystal=False,molecule=True)
           proc_adsorb(cryst,mol,data)

def cleave_operation():
    struct=readstructure()
    if isinstance(Structure,Molecule):
        print("cleave operation is only supported for periodic structure")
        return
    print('your choice ?')
    print('{} >>> {}'.format('1','cleave surface'))
    print('{} >>> {}'.format('2','cleave sphere cluster'))
    print('{} >>> {}'.format('3','cleave shell structure'))
    wait_sep()
    in_str=""
    while in_str=="":
       in_str=input().strip()
    choice=int(in_str)
    if choice==1:
        print(" input the miller index, minimum size in angstroms of layers containing atomssupercell")
        print(" and Minimize size in angstroms of layers containing vacuum like this:")
        print(' 1 0 0 | 5 | 5')
        print(' it means miller index is [1,0,0]')
        print(" min_slab_size is 5 Ang ")
        print(" min_vacum_size is 5 Ang ")
        print(" or like this : ")
        print(' 2 | 5 | 5')
        print(' it will generate all slab with miller index less than 2')
         
        def generate_selected_slab(in_str):
            tmp_list=in_str.split('|')
            miller_index=[int(x) for x in tmp_list[0].strip().split() ]
            min_slab_size=float(tmp_list[1])
            min_vac_size=float(tmp_list[2])
            slab=SlabGenerator(struct,miller_index,min_slab_size=min_slab_size,min_vacuum_size=min_vac_size,lll_reduce=True)
            slab_struct=slab.get_slab()
            slab_struct.sort()
            miller_str=[str(i) for i in miller_index]
            filename='_'.join(miller_str)+'.vasp'
            slab_struct.to(filename=filename,fmt='POSCAR')
        def generate_all_slab(in_str):
            tmp_list=in_str.split('|')
            max_index=int(tmp_list[0])
            min_slab_size=float(tmp_list[1])
            min_vac_size=float(tmp_list[2])
            slabs=generate_all_slabs(struct,max_index=max_index,min_slab_size=min_slab_size,min_vacuum_size=min_vac_size,lll_reduce=True)
            for slab_struct in slabs:
                slab_struct.sort()
                miller_str=[str(i) for i in slab_struct.miller_index]
                filename='_'.join(miller_str)+'.vasp'
                slab_struct.to(filename=filename,fmt='POSCAR')
        wait_sep()
        in_str=""
        while in_str=="":
           in_str=input().strip()
        len_para=len(in_str.split('|')[0].split())
        #if in_str.strip().startswith('a'):
        if len_para==3:
           generate_selected_slab(in_str)
        #elif in_str.strip().startswith('b'):
        elif len_para==1:
           generate_all_slab(in_str)
        else:
           print("unknow format")
           os._exit()
 
    elif choice==2:
        print(" input the center atom index, sphere radius and vacuum layer thickness")
        print(' 1 3.5 15')
        print(' it means the sphere will be selected according to the 1st atom')
        print(" with the radius equals 5Ang, and vacuum layer thickness is 15 Ang")
        wait_sep()
        in_str=""
        while in_str=="":
           in_str=input().strip()
        para=in_str.split()
        center_atom=int(para[0])-1 
        radius=float(para[1])
        vacuum=float(para[2])
        center_coord=struct[center_atom].coords
        sites=struct.get_neighbors_in_shell(center_coord,0,radius)
        coords=[site[0].coords for site in sites]
        species=[site[0].specie for site in sites]
        mol=Molecule(coords=coords,species=species)
        max_dist=np.max(mol.distance_matrix)
        a=b=c=max_dist+vacuum
        box_struct=mol.get_boxed_structure(a,b,c)
        file_name="sphere.vasp"
        box_struct.to(filename=file_name,fmt='poscar')
    elif choice==3:
        print(" input the center atom index, start radius, shell thickness and")
        print(" vacuum layer thickness")
        print(' 1 5 10  15')
        print(' it means the ball shell will be selected according to the 1st atom')
        print(" with the 5< r <15Ang, and vacuum layer thickness is 15 Ang")
        wait_sep()
        in_str=""
        while in_str=="":
           in_str=input().strip()
        para=in_str.split()
        center_atom=int(para[0])-1 
        radius=float(para[1])
        shell=float(para[2])
        vacuum=float(para[3])
        center_coord=struct[center_atom].coords
        sites=struct.get_neighbors_in_shell(center_coord,radius,shell)
        coords=[site[0].coords for site in sites]
        species=[site[0].specie for site in sites]
        mol=Molecule(coords=coords,species=species)
        max_dist=np.max(mol.distance_matrix)
        a=b=c=max_dist+vacuum
        box_struct=mol.get_boxed_structure(a,b,c)
        file_name="shell.vasp"
        box_struct.to(filename=file_name,fmt='poscar')

    else:
        print("unkown choice")
        return

def strain_operation():
    struct=readstructure()
    if isinstance(Structure,Molecule):
        print("strain operation is only supported for periodic structure")
        return
    print('input the strain component like :')
    print("0.01")
    print("it means aplly a strain of 1% along all directions")
    print("or")
    print("0.01 0.0 0.0")
    print("it means apply a strain of 1% along the x direction")
    print("or")
    print("0.01:0.03:5 0.0 0.0")
    print("it means to devide strain range into 5 parts")
    wait_sep()
    strain_str=input()
    tmp_list=strain_str.split(" ")
    if len(tmp_list)==1:
        strain=[[float(tmp_list[0]),float(tmp_list[0]),float(tmp_list[0])]]
    elif len(tmp_list)==3 and not ":" in strain_str:
        strain=[[float(x) for x in tmp_list]]
    elif len(tmp_list)==3 and strain_str.count(":")%2==0:
        tmp1=[float(x) for x in tmp_list[0].split(":")]
        tmp2=[float(x) for x in tmp_list[1].split(":")]
        tmp3=[float(x) for x in tmp_list[2].split(":")]
        strain=[]
        if len(tmp1)==3:
            x_range=np.linspace(tmp1[0],tmp1[1],int(tmp1[2]))
        else:
            x_range=tmp1

        if len(tmp2)==3:
            y_range=np.linspace(tmp2[0],tmp2[1],int(tmp2[2]))
        else:
            y_range=tmp2
        if len(tmp3)==3:
            z_range=np.linspace(tmp3[0],tmp3[1],int(tmp3[2]))
        else:
            z_range=tmp3
        for ix in x_range:
            for iy in y_range:
                for iz in z_range:
                        strain.append([ix,iy,iz])
    else:
        print("unknow format")
    i_count=0
    fis_name='index_strain.dat'
    fis=open(fis_name,'w')
    for i_strain in strain:
        outfile_name='strain_'+str(i_count)+'.vasp'
       # print(outfile_name)
        struct_cp=struct.copy()
        struct_cp.apply_strain(i_strain)
       # print(struct_cp)
        fis.writelines("%3d %7.4f %7.4f %7.4f\n"%(i_count,strain[i_count][0],strain[i_count][1],strain[i_count][2]))
        struct_cp.to(filename=outfile_name,fmt='poscar')
        i_count+=1
    fis.close()
#    print(strain)

def move_to_zcenter(struct):
    frac=struct.frac_coords
    zcenter=np.mean(frac[:,2])
    zshift=0.5-zcenter
    frac[:,2]=frac[:,2]+zshift
    return Structure(struct.lattice,struct.species,frac)
    
def twoD_operation():
    #sepline(ch='2D structure operation',sp='=')
    print('your choice ?')
    print('{} >>> {}'.format('1','build rippled structure'))
    print('{} >>> {}'.format('2','build multi-layered structure'))
    print('{} >>> {}'.format('3','split multi-layered structure'))
    print('{} >>> {}'.format('4','resize vacuum layer'))
    print('{} >>> {}'.format('5','center atomic-layer along z direction'))
    print('{} >>> {}'.format('6','apply strain along different direction'))
    print('{} >>> {}'.format('7','constrain atom in specific range'))
    print('{} >>> {}'.format('8','get a substrate for 2D material (online!!!)'))
    wait_sep()
    
    
    in_str=""
    while in_str=="":
          in_str=input().strip()
    
    choice=int(in_str)
    if choice==1:
        struct=readstructure()
        margin_dist=1.0 # 
        assert(struct.lattice.is_orthogonal)
        print("input the supercell scaling factor")
        print('for x direction can be: 10 1 1')
        print('for y direction can be: 1 10 1')
        wait_sep()
        scale_str=input()
        scaling=[int(x) for x in scale_str.split()]
        if (len(scaling)!=3):
            print('unknow format')
        else:
            if scaling[0]>=scaling[1]:
                direction=0  # for x direction
            else:
                direction=1  # for y direction
        print("input the strain range")
        print("example: 0.02:0.1:10 ")
        wait_sep()
        strain_str=input()
        tmp=[float(x) for x in strain_str.split(":")]
        strain=[]
        if len(tmp)==3:
             strain_range=np.linspace(tmp[0],tmp[1],int(tmp[2]))
        #     print(strain_range)
        else:
            print("unknow format")
            return
        print("input the index of atom need to be fixed")
        print("example: 1 10 11 20 ")
        print("0 means fix the atom automatically")
        wait_sep()
        atom_index_str=input()
        try:
            atom_index=[int(x) for x in strain_str.split("")]
            auto_fix=False
        except:
            atom_index=[int(atom_index_str)]
            auto_fix=True
        
        struct_sc=struct.copy()
        struct_sc.make_supercell(scaling)
        natom=struct_sc.num_sites
        min_z=np.min(struct_sc.cart_coords[:,2])
        tmp_coords=np.ones((struct_sc.num_sites,1))*min_z
        cart_coords=struct_sc.cart_coords;
        cart_coords[:,2]=cart_coords[:,2]-tmp_coords.T+0.01
        frac_coords_new=np.dot(cart_coords,np.linalg.inv(struct_sc.lattice.matrix))
        for i_strain in strain_range:
            new_lat_matrix=struct_sc.lattice.matrix
            new_lat_matrix[direction,direction]=struct_sc.lattice.matrix[direction,direction]*(1-i_strain)
            fname="%10.5f"%(i_strain)+'_wo.vasp' # structure only applied with in-plan strain
            struct_wo_ripple=Structure(new_lat_matrix,struct_sc.species,frac_coords_new)
            struct_wo_ripple.to(filename=fname.strip(),fmt='poscar')
            frac_coords_new_cp=struct_wo_ripple.frac_coords.copy()
            cart_coords_new_cp=struct_wo_ripple.cart_coords.copy()
            nz=0
            selective_dynamics=[[True for col in range(3)] for row in range(natom)]  
            z_shift=np.zeros((natom,3))
            for i_atom in range(natom):
                z_shift[i_atom,2]=40*(2*i_strain-10*i_strain**2)*np.sin(cart_coords_new_cp[i_atom,direction]*np.pi/new_lat_matrix[direction,direction]) 
                if cart_coords_new_cp[i_atom,direction]< nz or cart_coords_new_cp[i_atom,direction] > new_lat_matrix[direction,direction] -nz:
                   z_shift[i_atom,2]=0.0
                   
                if auto_fix:
                   if struct_wo_ripple[i_atom].coords[direction]<margin_dist or \
                      struct_wo_ripple[i_atom].coords[direction] > new_lat_matrix[direction,direction]-margin_dist:
                      selective_dynamics[i_atom]=[False,False,False]
                else:     
                   if i_atom in atom_index:
                      selective_dynamics[i_atom]=[False,False,False]
            
            struct_w_ripple=Structure(new_lat_matrix,struct_sc.species,cart_coords_new_cp+z_shift,coords_are_cartesian=True,\
                                     site_properties={'selective_dynamics':selective_dynamics})
            fname="%10.5f"%(i_strain)+'_w.vasp' # structure  applied with in-plan strain  and ripple    
            struct_w_ripple.to(filename=fname.strip(),fmt='poscar')        
    elif choice==2:
        struct=readstructure()
        print("input the number of layers")
        wait_sep()
        in_str=""
        while in_str=="":
              in_str=input().strip()
        layer_number=int(in_str)

        print("input the layer distance")
        wait_sep()
        species=[]
        in_str=""
        while in_str=="":
              in_str=input().strip()
        layer_distance=float(in_str)

        new_struct=move_to_zcenter(struct)
        struct_thickness=np.max(new_struct.cart_coords[:,2])-np.min(new_struct.cart_coords[:,2])
        natom=new_struct.num_sites
        new_cart_coords=np.zeros((natom*layer_number,3))
        for i in range(layer_number):
            new_cart_coords[i*natom:(i+1)*natom,0:2]=new_struct.cart_coords[:,0:2]
            new_cart_coords[i*natom:(i+1)*natom,2]=new_struct.cart_coords[:,2]+i*(layer_distance+struct_thickness)
            species.extend(new_struct.species)
        new_lat=new_struct.lattice.matrix
        new_lat[2,2]=new_lat[2,2]+layer_distance*layer_number
        tmp_struct=Structure(new_lat,species,new_cart_coords,coords_are_cartesian=True)
        tmp1_struct=move_to_zcenter(tmp_struct)
        tmp2_struct=tmp1_struct.get_sorted_structure()
        tmp2_struct.to(filename='layer_'+str(layer_number)+'.vasp',fmt='poscar')
 
    elif choice==3:
        ProperDist=3.5  # Ang
        struct=readstructure()
        (atom_index,in_str)=atom_selection(struct)
        print("input the splitting distance, 10 Ang is enough!")
        wait_sep()
        in_str=""
        while in_str=="":
              in_str=input().strip()
        SplitDistance=float(in_str)
        print("numbers of splitting site, 50 sites are enough!")
        wait_sep()
        in_str=""
        while in_str=="":
              in_str=input().strip()
        NumberSplitSite=int(in_str)
        DensityN=int(NumberSplitSite*0.75)
        SparseN=NumberSplitSite-DensityN+1
#        print(DensityN,SparseN)
        dist=ProperDist/(DensityN-1)
        SplitDistanceArray=np.zeros(NumberSplitSite+1)        
        for Nsite in range(DensityN):
            SplitDistanceArray[Nsite]=(Nsite)*dist
        
        dist=(SplitDistance-ProperDist)/SparseN
        for  Nsite in range(SparseN):
             SplitDistanceArray[Nsite+DensityN]=ProperDist+(Nsite+1)*dist

#        print(SplitDistanceArray)
        coords=struct.cart_coords
        for Nsite in range(NumberSplitSite+1):
            coords=struct.cart_coords
            for atom in atom_index:
                coords[atom,2]=coords[atom,2]+SplitDistanceArray[Nsite]
            tmp_struct=Structure(struct.lattice,struct.species,coords,coords_are_cartesian=True) 
            fname=str(Nsite)+'.vasp'
            tmp_struct.to(filename=fname,fmt='poscar')
        data=np.zeros((NumberSplitSite+1,2))
        for i,j in enumerate(SplitDistanceArray):
            data[i][0]=i
            data[i][1]=j
        head_line="#%(key1)+12s  %(key2)+12s"%{'key1':'index','key2':'distance/Ang'}
        fmt="%12d %12.6f"+'\n'
        write_col_data('split.dat',data,head_line,sp_fmt=fmt) 
        return

    elif choice==4:
        struct=readstructure()
        new_struct=move_to_zcenter(struct)
        struct_thickness=np.max(new_struct.cart_coords[:,2])-np.min(new_struct.cart_coords[:,2])
        vac_layer_thickness=new_struct.lattice.c-struct_thickness
        print("current vacuum layer thickness is %6.3f Ang"%(vac_layer_thickness))
        print("input the new value of vacuum layer thickness")
        wait_sep()
        in_str=""
        while in_str=="":
              in_str=input().strip()
        nvac_layer_thickness=float(in_str)

        assert(nvac_layer_thickness>struct_thickness)
        new_lat=new_struct.lattice.matrix
        new_lat[2,2]=nvac_layer_thickness
        tmp_struct=Structure(new_lat,new_struct.species,new_struct.cart_coords,coords_are_cartesian=True)
        center_struct=move_to_zcenter(tmp_struct)
        center_struct.to(filename='new_vacuum.vasp',fmt='poscar') 
        return

    elif choice==5:
        struct=readstructure()
        new_struct=move_to_zcenter(struct)
        new_struct.to(filename='z-center.vasp',fmt='poscar')
        return

    elif choice==6:
        struct=readstructure()
        assert(struct.lattice.is_orthogonal)
        try:
           import sympy
        except:
           print("you must install sympy module")
           return

        new_struct=move_to_zcenter(struct)
        print("input the elastic of material by order : C11 C12 C22 C66")
        wait_sep()
        in_str=""
        while in_str=="":
              in_str=input().strip().split()
        elastic_constant=[float(x) for x in in_str] 
        if len(elastic_constant) !=4:
           print("you must input C11 C12 C22 C66")
           return
        C11=elastic_constant[0]
        C12=elastic_constant[1]
        C22=elastic_constant[2]
        C66=elastic_constant[3]
        print("input applied force: e.x. 1.0 GPa nm")
        wait_sep()
        in_str=""
        while in_str=="":
              in_str=input().strip()
        sigma=float(in_str)
        orig_struct=new_struct.copy()
        new_struct =new_struct.copy()
        natom=orig_struct.num_sites
        lat=orig_struct.lattice.matrix
        pos=orig_struct.frac_coords
        nps=37
        phi=np.linspace(0,360,nps)*np.pi/180
        vzz=C12/C22
        temp_num=(C11*C22-C12**2)/(C22*C66)
        d1=C11/C22 +1.0-temp_num;
        d2=-(2.0*C12/C22-temp_num);
        d3=C11/C22;
        F=sigma*C22/(C11*C22-C12**2.0); 
        Poisson=(vzz*(np.cos(phi))**4.0-d1*(np.cos(phi))**2.0*(np.sin(phi))**2.0+vzz*(np.sin(phi))**4.0)/\
                ((np.cos(phi))**4.0+d2*(np.cos(phi))**2.0*(np.sin(phi))**2.0+d3*(np.sin(phi))**4.0)
        
        eps_theta=F*((np.cos(phi))**4+d2*(np.cos(phi))**2.0*(np.sin(phi))**2.0+d3*(np.sin(phi))**4.0) 
        t = sympy.Symbol('t', real=True)
        e = sympy.Symbol('e', real=True)
        v = sympy.Symbol('v', real=True)
        eprim=sympy.Matrix([[e+1,0],[ 0,1-e*v]])
        R=sympy.Matrix([[sympy.cos(t),-sympy.sin(t)],[ sympy.sin(t),sympy.cos(t)]])
        eps_mat=R*eprim*R.adjugate()
        for k in range(len(phi)):
            cur__phi=phi[k]*180/np.pi 
            Rot=eps_mat.subs({e:eps_theta[k],v:Poisson[k],t:phi[k]})
            fname=str(k)+'.vasp'
            final_lat=np.matrix(np.eye(3))
            final_lat[0,0]=Rot[0,0]
            final_lat[0,1]=Rot[0,1]
            final_lat[1,0]=Rot[1,0]
            final_lat[1,1]=Rot[1,1]
            lat_new=lat*final_lat
            tmp_struct=Structure(lat_new,new_struct.species,pos)
            tmp_struct.to(filename=fname,fmt='poscar') 
        return

    elif choice==7:
        struct=readstructure()
        natom=struct.num_sites
        atom_index,in_str=atom_selection(struct)
        selective_dynamics=[[True for col in range(3)] for row in range(natom)]  
        for i in range(natom):
               if i in atom_index:
                    selective_dynamics[i]=[False,False,False]
        tmp_struct=Structure(struct.lattice,struct.species,struct.frac_coords,site_properties={'selective_dynamics':selective_dynamics})
        poscar=Poscar(tmp_struct)
        poscar.comment=poscar.comment+' |--> '+in_str
        poscar.write_file('Fixed.vasp')
        return

    elif choice==8:
        print('your choice ?')
        print('{} >>> {}'.format('1','input 2D structure from local disk'))
        print('{} >>> {}'.format('2','get 2D structure online'))
        wait_sep()
        in_str=""
        while in_str=="":
              in_str=input().strip()
        choice=int(in_str)

        if choice==1:
           mpid=None
           struct=readstructure()
        else:
           print("input the mp-id for your structure")
           wait_sep()
           in_str=""
           while in_str=="":
                 in_str=input().strip()
           mpid=in_str
           struct=None
          
        film,substrates=make_connect(mpid=mpid,struct=struct)
        df=get_subs(film,substrates)
        df.to_csv('substrate.csv', sep=',', header=True, index=True)
        return
    else:
        print("unkonw choice")
        return

def get_xrd():
    struct=readstructure()
    xrd=XRDCalculator()
  
#    xrd_data=xrd.get_xrd_pattern(struct)
    xrd_data=xrd.get_pattern(struct)
    jxrd_data=jsanitize(xrd_data.as_dict())
    fname='XRD.json'
    proc_str="Saving data to "+ fname +" File ..."
    procs(proc_str,1,sp='-->>')
    json_store(jxrd_data,fname)

    data=np.vstack((xrd_data.x,xrd_data.y)).T
    
    margin=10.
    ds=xrd_data.x[0]-margin
    de=xrd_data.x[-1]+margin
    tmp_data=[ds]+xrd_data.x.tolist()+[de]
    tmp_data1=np.diff(tmp_data).tolist()
    tmp_data2=np.array([0]+np.cumsum(tmp_data1).tolist())
    tmp_data3=tmp_data2/tmp_data2[-1]
    x_data=np.linspace(ds,de,10000)
    y_data=np.zeros((len(x_data)))
    for i in range(1,len(tmp_data3)-1):
      index=int(tmp_data3[i]*10000)
      y_data[index]=xrd_data.y[i-1]  

    data=np.vstack((x_data,y_data))
    data=(smear(data, sigma=0.1)).T

    head_line="#%(key1)+12s %(key2)+12s"%{'key1':'2theta','key2':'Intensity'}
    fname='XRD.dat'
    proc_str="Saving data to "+ fname +" File ..."
    procs(proc_str,2,sp='-->>')
    write_col_data(fname,data,head_line)

    check_matplotlib()
    fname='XRD.png'
    proc_str="Saving plot to "+ fname +" File ..."
    procs(proc_str,3,sp='-->>')
    plt=xrd.get_plot(struct)
    plt.savefig(fname,format='png')

