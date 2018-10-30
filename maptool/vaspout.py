#!/usr/bin/env python
from    __future__ import division, unicode_literals, print_function
import os
import sys
import itertools
import numpy as np
from   numpy import linalg as LA
from   copy import copy
from   pymatgen import Structure
from   pymatgen.io.vasp import VolumetricData,Poscar,Vasprun,Outcar,Locpot,Procar
from   pymatgen.electronic_structure.core import Spin
from   monty.io import zopen, reverse_readfile
from   pymatgen.electronic_structure.bandstructure import BandStructure, \
       BandStructureSymmLine, get_reconstructed_band_structure
from   pymatgen.electronic_structure.core import Spin, Orbital, OrbitalType, Magmom
from   pymatgen.electronic_structure.dos import CompleteDos, Dos
from   pymatgen.electronic_structure.plotter import DosPlotter,BSPlotter
from   pymatgen.symmetry import analyzer
from   maptool.constants import H,C0,Avogadro
from   maptool.proc_mech import calc_elastic_prop
from   maptool.data_io import write_col_data
from   maptool.utils import atom_selection,procs,wait_sep,check_file,check_matplotlib

min_gap=0.001
min_mag=0.001

def total_dos():
   check_matplotlib()
   filename='vasprun.xml'
   check_file(filename)
   proc_str="Reading Data From "+ filename +" File ..."
   procs(proc_str,1,sp='-->>')
   vsr=Vasprun(filename,parse_eigen=False)
   tdos=vsr.tdos
   idos=vsr.idos
   E=tdos.energies-tdos.efermi
   if vsr.is_spin:
      proc_str="This Is a Spin-polarized Calculation."
      procs(proc_str,0,sp='-->>')
      proc_str="Writting TDOS.dat File ..."
      TDOSUP=tdos.densities[Spin.up]
      TDOSDOWN=tdos.densities[Spin.down]
      ETDOS=np.vstack((E,TDOSUP,TDOSDOWN))
      head_line="#%(key1)+12s%(key2)+12s%(key3)+12s"%{'key1':'Energy(eV)','key2':'SpinUp','key3':'SpinDown'}
      write_col_data('TDOS.dat',ETDOS.T,head_line)

      proc_str="Writting IDOS.dat File ..."
      procs(proc_str,3,sp='-->>')
      IDOSUP=idos.densities[Spin.up]
      IDOSDOWN=idos.densities[Spin.down]
      EIDOS=np.vstack((E,IDOSUP,IDOSDOWN))
      head_line="#%(key1)+12s%(key2)+12s%(key3)+12s"%{'key1':'Energy(eV)','key2':'IntSpinUp','key3':'IntSpinDown'}
      write_col_data('IDOS.dat',EIDOS.T,head_line)

      plt1=DosPlotter()
      plt2=DosPlotter()
      plt1.add_dos('Total DOS',tdos)
      plt2.add_dos('Total DOS',idos)
      try:
       # plt1.show()
        plt1.save_plot('TotalDOS.png', img_format="png")
       # plt2.show()
        plt2.save_plot('IntegratedDOS.png', img_format="png")
      except:
        print("pls use gnuplot to plot TDOS.dat and IDOS.dat")
   else:
      if vsr.parameters['LNONCOLLINEAR']:
         proc_str="This Is a Non-Collinear Calculation."
      else:
          proc_str="This Is a Non-Spin Calculation."
      procs(proc_str,0,sp='-->>')
      proc_str="Writting TDOS.dat File ..."
      procs(proc_str,2,sp='-->>')
      
      TDOS=tdos.densities[Spin.up]
      ETDOS=np.vstack((E,TDOS))
      head_line="#%(key1)+12s%(key2)+12s"%{'key1':'Energy(eV)','key2':'TotalDOS'}
      write_col_data('TDOS.dat',ETDOS.T,head_line)

      proc_str="Writting IDOS.dat File ..."
      procs(proc_str,3,sp='-->>')
      IDOS=idos.densities[Spin.up]
      EIDOS=np.vstack((E,IDOS))
      head_line="#%(key1)+12s%(key2)+13s"%{'key1':'Energy(eV)','key2':'IntegratedDOS'}
      write_col_data('IDOS.dat',EIDOS.T,head_line)

      plt1=DosPlotter()
      plt2=DosPlotter()
      plt1.add_dos('Total DOS',tdos)
      plt2.add_dos('Integrated DOS',idos)
      filename4="TotalDOS.png IntegratedDOS.png"
      proc_str="Saving Plot to "+ filename4 +" File ..."
      procs(proc_str,4,sp='-->>')

      try:
       # plt1.show()
        plt1.save_plot('TotalDOS.png', img_format="png")
       # plt2.show()
        plt2.save_plot('IntegratedDOS.png', img_format="png")
      except:
        print("pls use gnuplot to plot TDOS.dat and IDOS.dat")
  
def projected_dos():
   step_count=1
   filename='vasprun.xml'
   check_file(filename)
   proc_str="Reading Data From "+ filename +" File ..."
   procs(proc_str,step_count,sp='-->>')
   vsr=Vasprun(filename)
   nedos=vsr.parameters['NEDOS']
   struct=vsr.final_structure
   pdos=vsr.pdos

   filename='PROCAR'
   check_file(filename)
   step_count+=1
   proc_str="Reading Data From "+ filename +" File ..."
   procs(proc_str,step_count,sp='-->>')
   procar=Procar(filename)
   nbands=procar.nbands
   nions=procar.nions
   norbitals=len(procar.orbitals)
   nkpoints=procar.nkpoints

   (atom_index,in_str)=atom_selection(struct)

   if len(atom_index)==0:
      print("No atoms selected!")
      return
#   print(atom_index)

   if vsr.is_spin:
      proc_str="This Is a Spin-polarized Calculation."
      procs(proc_str,0,sp='-->>')

      contrib=np.zeros((nedos,norbitals+1,2))
      energies=vsr.tdos.energies-vsr.efermi
      for ispin in [0,1]:
          if ispin==0:
              spin=Spin.up
              s_name='Up'
          else:
              spin=Spin.down
              s_name='Down'

          contrib[:,0,ispin]=energies  
          for i in atom_index:
                for j in range(norbitals):
                     contrib[:,j+1,ispin]=contrib[:,j+1,ispin]+pdos[i][Orbital(j)][spin]

          step_count+=1
          filename="PDOS_"+s_name+".dat"
          proc_str="Writting Projected DOS Data to "+ filename +" File ..."
          procs(proc_str,step_count,sp='-->>')
          tmp1_str="#%%(key1)+12s"
          tmp2_dic={'key1':'Energy(ev)'}
          for i in range(norbitals):
              tmp1_str+="%(key"+str(i+2)+")+12s"
              tmp2_dic["key"+str(i+2)]=procar.orbitals[i]

#          print(tmp1_str)
          atom_index_str=[str(x+1) for x in atom_index]
          head_line1="#String: "+in_str+'\n#Selected atom: ' +' '.join(atom_index_str)+'\n'
          head_line2=tmp1_str % tmp2_dic
          head_line=head_line1+head_line2
          write_col_data(filename,contrib[:,:,ispin],head_line)

   else:
      if vsr.parameters['LNONCOLLINEAR']:
         proc_str="This Is a Non-Collinear Calculation."
         procs(proc_str,0,sp='-->>')
      else:
         proc_str="This Is a Non-Spin Calculation."
         procs(proc_str,0,sp='-->>')

      contrib=np.zeros((nedos,norbitals+1))
      energies=vsr.tdos.energies-vsr.efermi
      contrib[:,0]=energies
      for i in atom_index:
            for j in range(norbitals):
                 contrib[:,j+1]=contrib[:,j+1]+pdos[i][Orbital(j)][Spin.up]

      step_count+=1
      filename="PDOS.dat"
      proc_str="Writting Projected DOS Data to "+ filename +" File ..."
      procs(proc_str,step_count,sp='-->>')
      tmp1_str="#%(key1)+12s%(key2)+12s"
      tmp2_dic={'key1':'K-Distance','key2':'Energy(ev)'}
      for i in range(norbitals):
          tmp1_str+="%(key"+str(i+3)+")+12s"
          tmp2_dic["key"+str(i+3)]=procar.orbitals[i]

#      print(tmp1_str)
      atom_index_str=[str(x+1) for x in atom_index]
      head_line1="#String: "+in_str+'\n#Selected atom: ' +' '.join(atom_index_str)+'\n'
      head_line2=tmp1_str % tmp2_dic
      head_line=head_line1+head_line2
      write_col_data(filename,contrib,head_line)

def band_structure():
    check_matplotlib()
    step_count=1

    filename='vasprun.xml'
    check_file(filename)
    proc_str="Reading Data From "+ filename +" File ..."
    procs(proc_str,step_count,sp='-->>')
    vsr=Vasprun(filename)

    step_count+=1
    filename='KPOINTS'
    check_file(filename)
    proc_str="Reading Data From "+ filename +" File ..."
    procs(proc_str,step_count,sp='-->>')
    bands = vsr.get_band_structure(filename, line_mode=True, efermi=vsr.efermi)

    step_count+=1
    filename='OUTCAR'
    check_file(filename)
    proc_str="Reading Data From "+ filename +" File ..."
    procs(proc_str,step_count,sp='-->>')
    outcar=Outcar('OUTCAR')
    mag=outcar.as_dict()['total_magnetization']

    if vsr.is_spin:
       proc_str="This Is a Spin-polarized Calculation."
       procs(proc_str,0,sp='-->>')
       tdos=vsr.tdos
       SpinUp_gap=tdos.get_gap(spin=Spin.up) 
       cbm_vbm_up=tdos.get_cbm_vbm(spin=Spin.up)
       SpinDown_gap=tdos.get_gap(spin=Spin.down) 
       cbm_vbm_down=tdos.get_cbm_vbm(spin=Spin.up)

       if SpinUp_gap > min_gap and SpinDown_gap > min_gap:
          is_metal=False
          is_semimetal=False
       elif SpinUp_gap > min_gap and SpinDown_gap < min_gap:
          is_metal=False
          is_semimetal=True
       elif SpinUp_gap < min_gap and SpinDown_gap > min_gap:
          is_metal=False
          is_semimetal=True
       elif SpinUp_gap < min_gap and SpinDown_gap < min_gap:
          is_metal=True
          is_semimetal=False
          
       if is_metal:   
          proc_str="This Material Is a Metal."
          procs(proc_str,0,sp='-->>')
       if not is_metal and is_semimetal:
          proc_str="This Material Is a Semimetal."
          procs(proc_str,0,sp='-->>')
       else:
          proc_str="This Material Is a Semiconductor."
          procs(proc_str,0,sp='-->>')
          proc_str="Total magnetization is "+str(mag)
          procs(proc_str,0,sp='-->>')
          if mag > min_mag:
             proc_str="SpinUp  : vbm=%f eV cbm=%f eV gap=%f eV"%(cbm_vbm_up[1],cbm_vbm_up[0],SpinUp_gap)
             procs(proc_str,0,sp='-->>')
             proc_str="SpinDown: vbm=%f eV cbm=%f eV gap=%f eV"%(cbm_vbm_down[1],cbm_vbm_down[0],SpinUp_gap)
             procs(proc_str,0,sp='-->>')
          else:
             proc_str="SpinUp  : vbm=%f eV cbm=%f eV gap=%f eV"%(cbm_vbm_up[1],cbm_vbm_up[0],SpinUp_gap)
             procs(proc_str,0,sp='-->>')
       step_count+=1
       filename="BAND.dat"
       proc_str="Writting Band Structure Data to "+ filename +" File ..."
       procs(proc_str,step_count,sp='-->>')
       band_data_up=bands.bands[Spin.up]
       band_data_down=bands.bands[Spin.down]
       y_data_up=band_data_up.reshape(1,band_data_up.shape[0]*band_data_up.shape[1])[0]-vsr.efermi #shift fermi level to 0
       y_data_down=band_data_down.reshape(1,band_data_down.shape[0]*band_data_down.shape[1])[0]-vsr.efermi #shift fermi level to 0
       x_data=np.array(bands.distance*band_data_up.shape[0])
       data=np.vstack((x_data,y_data_up,y_data_down)).T
       head_line="#%(key1)+12s%(key2)+13s%(key3)+15s"%{'key1':'K-Distance','key2':'UpEnergy(ev)','key3':'DownEnergy(ev)'}
       write_col_data(filename,data,head_line,band_data_up.shape[1])
 
    else:
       if vsr.parameters['LNONCOLLINEAR']:
          proc_str="This Is a Non-Collinear Calculation."
       else:
           proc_str="This Is a Non-Spin Calculation."
       procs(proc_str,0,sp='-->>')
       cbm=bands.get_cbm()['energy']
       vbm=bands.get_vbm()['energy']
       gap=bands.get_band_gap()['energy']
       if not bands.is_metal():
          proc_str="This Material Is a Semiconductor."
          procs(proc_str,0,sp='-->>')
          proc_str="vbm=%f eV cbm=%f eV gap=%f eV"%(vbm,cbm,gap)
          procs(proc_str,0,sp='-->>')
       else:
          proc_str="This Material Is a Metal."
          procs(proc_str,0,sp='-->>')
       
       step_count+=1
       filename3="BAND.dat"
       proc_str="Writting Band Structure Data to "+ filename3 +" File ..."
       procs(proc_str,step_count,sp='-->>')
       band_data=bands.bands[Spin.up]
       y_data=band_data.reshape(1,band_data.shape[0]*band_data.shape[1])[0]-vsr.efermi #shift fermi level to 0
       x_data=np.array(bands.distance*band_data.shape[0])
       data=np.vstack((x_data,y_data)).T
       head_line="#%(key1)+12s%(key2)+13s"%{'key1':'K-Distance','key2':'Energy(ev)'}
       write_col_data(filename3,data,head_line,band_data.shape[1])
       step_count+=1
       bsp=BSPlotter(bands)
       filename4="HighSymmetricPoints.dat"
       proc_str="Writting Label infomation to "+ filename4 +" File ..."
       procs(proc_str,step_count,sp='-->>')
       head_line="#%(key1)+12s%(key2)+12s%(key3)+12s"%{'key1':'index','key2':'label','key3':'position'}
       line=head_line+'\n'
       for i,label in enumerate(bsp.get_ticks()['label']):
           new_line="%(key1)12d%(key2)+12s%(key3)12f\n"%{'key1':i,'key2':label,'key3':bsp.get_ticks()['distance'][i]}
           line+=new_line
       line+='\n'
       write_col_data(filename4,line,'',str_data=True) 
    try:
       step_count+=1
       filename5="BAND.png"
       proc_str="Saving Plot to "+ filename5 +" File ..."
       procs(proc_str,step_count,sp='-->>')
       bsp.save_plot(filename5, img_format="png")
    except:
       print("Figure output fails !!!")   

      
def projected_band_structure():
   step_count=1
   filename='vasprun.xml'
   check_file(filename)
   proc_str="Reading Data From "+ filename +" File ..."
   procs(proc_str,step_count,sp='-->>')
   vsr=Vasprun(filename)

   filename='PROCAR'
   check_file(filename)
   step_count+=1
   proc_str="Reading Data From "+ filename +" File ..."
   procs(proc_str,step_count,sp='-->>')
   procar=Procar(filename)
   nbands=procar.nbands
   nions=procar.nions
   norbitals=len(procar.orbitals)
   nkpoints=procar.nkpoints

   step_count+=1
   filename='KPOINTS'
   check_file(filename)
   proc_str="Reading Data From "+ filename +" File ..."
   procs(proc_str,step_count,sp='-->>')
   bands = vsr.get_band_structure(filename, line_mode=True, efermi=vsr.efermi)
   struct=vsr.final_structure
   (atom_index,in_str)=atom_selection(struct)
   
   if len(atom_index)==0:
      print("No atoms selected!")
      return
#   print(atom_index)

   if vsr.is_spin:
      proc_str="This Is a Spin-polarized Calculation."
      procs(proc_str,0,sp='-->>')
      ISPIN=2
      contrib=np.zeros((nkpoints,nbands,norbitals,2))
      for i in atom_index:
          contrib[:,:,:,0]=contrib[:,:,:,0]+procar.data[Spin.up][:,:,i,:]
          contrib[:,:,:,1]=contrib[:,:,:,1]+procar.data[Spin.down][:,:,i,:]

      for ispin in range(2):
          proj_band=contrib[:,:,:,ispin].reshape(nkpoints*nbands,norbitals)
          step_count+=1
          if ispin==0:
              filename="PBAND_Up.dat"
          else:
              filename="PBAND_Down.dat"
          proc_str="Writting Projected Band Structure Data to "+ filename +" File ..."
          procs(proc_str,step_count,sp='-->>')
          band_data=bands.bands[Spin.up]
          y_data=band_data.reshape(1,nbands*nkpoints)[0]-vsr.efermi #shift fermi level to 0
          x_data=np.array(bands.distance*nbands)
          data=np.vstack((x_data,y_data,proj_band.T)).T
          tmp1_str="#%(key1)+12s%(key2)+12s"
          tmp2_dic={'key1':'K-Distance','key2':'Energy(ev)'}
          for i in range(norbitals):
              tmp1_str+="%(key"+str(i+3)+")+12s"
              tmp2_dic["key"+str(i+3)]=procar.orbitals[i]

#          print(tmp1_str)
          atom_index_str=[str(x+1) for x in atom_index]
          head_line1="#String: "+in_str+'\n#Selected atom: ' +' '.join(atom_index_str)+'\n'
          head_line2=tmp1_str % tmp2_dic
          head_line=head_line1+head_line2
          write_col_data(filename,data,head_line,nkpoints)

   else:
      if vsr.parameters['LNONCOLLINEAR']:
         proc_str="This Is a Non-Collinear Calculation."
         procs(proc_str,0,sp='-->>')
         ISPIN=3
      else:
         proc_str="This Is a Non-Spin Calculation."
         procs(proc_str,0,sp='-->>')
         ISPIN=1

      contrib=np.zeros((nkpoints,nbands,norbitals))
      for i in atom_index:
          contrib[:,:,:]=contrib[:,:,:]+procar.data[Spin.up][:,:,i,:]
 
      proj_band=contrib.reshape(nkpoints*nbands,norbitals)        
      step_count+=1
      filename="PBAND.dat"
      proc_str="Writting Projected Band Structure Data to "+ filename +" File ..."
      procs(proc_str,step_count,sp='-->>')
      band_data=bands.bands[Spin.up]
      y_data=band_data.reshape(1,nbands*nkpoints)[0]-vsr.efermi #shift fermi level to 0
      x_data=np.array(bands.distance*nbands)
      data=np.vstack((x_data,y_data,proj_band.T)).T
      tmp1_str="#%(key1)+12s%(key2)+12s"
      tmp2_dic={'key1':'K-Distance','key2':'Energy(ev)'}
      for i in range(norbitals):
          tmp1_str+="%(key"+str(i+3)+")+12s"
          tmp2_dic["key"+str(i+3)]=procar.orbitals[i]

#      print(tmp1_str)
      atom_index_str=[str(x+1) for x in atom_index]
      head_line1="#String: "+in_str+'\n#Selected atom: ' +' '.join(atom_index_str)+'\n'
      head_line2=tmp1_str % tmp2_dic
      head_line=head_line1+head_line2
      write_col_data(filename,data,head_line,nkpoints)

   step_count+=1
   bsp=BSPlotter(bands)
   filename="HighSymmetricPoints.dat"
   proc_str="Writting Label infomation to "+ filename +" File ..."
   procs(proc_str,step_count,sp='-->>')
   head_line="#%(key1)+12s%(key2)+12s%(key3)+12s"%{'key1':'index','key2':'label','key3':'position'}
   line=head_line+'\n'
   for i,label in enumerate(bsp.get_ticks()['label']):
       new_line="%(key1)12d%(key2)+12s%(key3)12f\n"%{'key1':i,'key2':label,'key3':bsp.get_ticks()['distance'][i]}
       line+=new_line
   write_col_data(filename,line,'',str_data=True) 
         
def select_one_band_structure():
   check_matplotlib()
   step_count=1
   filename='vasprun.xml'
   check_file(filename)
   proc_str="Reading Data From "+ filename +" File ..."
   procs(proc_str,step_count,sp='-->>')
   vsr=Vasprun(filename)

   step_count+=1
   filename='KPOINTS'
   check_file(filename)
   proc_str="Reading Data From "+ filename +" File ..."
   procs(proc_str,step_count,sp='-->>')
   bands = vsr.get_band_structure(filename, line_mode=True, efermi=vsr.efermi)
   nelect=vsr.parameters['NELECT']
   nbands=bands.nb_bands
   if vsr.is_spin:
      proc_str="This Is a Spin-polarized Calculation."
      procs(proc_str,0,sp='-->>')
      ISPIN=2
   else:
      if vsr.parameters['LNONCOLLINEAR']:
         proc_str="This Is a Non-Collinear Calculation."
         procs(proc_str,0,sp='-->>')
         ISPIN=3
      else:
         proc_str="This Is a Non-Spin Calculation."
         procs(proc_str,0,sp='-->>')
         ISPIN=1
   proc_str="Total band number is "+str(nbands)
   procs(proc_str,0,sp='-->>')
   proc_str="Total electron number is "+str(nelect)
   procs(proc_str,0,sp='-->>')

   print("which band would like to select ?")
   wait_sep()
   in_str=""
   while in_str=="":
         in_str=input().strip()
   selected_band=int(in_str)

   step_count+=1
   filename="BAND_"+str(selected_band)+'.dat'
   proc_str="Writting Selected Band Structure Data to "+ filename +" File ..."
   procs(proc_str,step_count,sp='-->>')
   if ISPIN==1 or ISPIN==3:
      band_data=bands.bands[Spin.up][selected_band-1]-vsr.efermi
      data=np.vstack((bands.distance,band_data)).T
      head_line="#%(key1)+12s%(key2)+13s"%{'key1':'K-Distance','key2':'Energy(ev)'}
      write_col_data(filename,data,head_line,len(band_data))
   else:
      band_data_up=bands.bands[Spin.up][selected_band-1]-vsr.efermi
      band_data_down=bands.bands[Spin.down][selected_band-1]-vsr.efermi
      data=np.vstack((bands.distance,band_data_up,band_data_down)).T
      head_line="#%(key1)+12s%(key2)+13s%(key3)+15s"%{'key1':'K-Distance','key2':'UpEnergy(ev)','key3':'DownEnergy(ev)'}
      write_col_data(filename,data,head_line,len(band_data_up))
   return

def charge_density():
   proc_str="Reading Charge Density From CHG File ..."
   procs(proc_str,1,sp='-->>')
   chg = CHGCAR.from_file("CHG")
   proc_str="Writing Charge Density To CHARGE.vasp File ..."
   procs(proc_str,2,sp='-->>')
   chg.write_file('CHGARGE.vasp','total')

class CHGCAR(VolumetricData):
    """
    Simple object for reading a CHGCAR file.

    Args:
        poscar (Poscar): Poscar object containing structure.
        data: Actual data.
    """

    def __init__(self, poscar, data, data_aug=[]):
        super(CHGCAR, self).__init__(poscar.structure, data, data_aug=data_aug)
        self.poscar = poscar
        self.name = poscar.comment
        self._distance_matrix = {}

    @staticmethod
    def from_file(filename):
        (poscar, data, data_aug) = VolumetricData.parse_file(filename)
        return CHGCAR(poscar, data, data_aug=data_aug)

    @property
    def net_magnetization(self):
        if self.is_spin_polarized:
            return np.sum(self.data['diff'])
        else:
            return None
    def write_file(self,file_name,data_type,vasp4_compatible=False):
        """
        Write the VolumetricData object to a vasp compatible file.

        Args:
            file_name (str): Path to a file
            vasp4_compatible (bool): True if the format is vasp4 compatible
        """

        def _print_fortran_float(f):
            """
            Fortran codes print floats with a leading zero in scientific
            notation. When writing CHGCAR files, we adopt this convention
            to ensure written CHGCAR files are byte-to-byte identical to
            their input files as far as possible.
            :param f: float
            :return: str
            """
            s = "{:.10E}".format(f)
            if f > 0:
                return "0."+s[0]+s[2:12]+'E'+"{:+03}".format(int(s[13:])+1)
            else:
                return "-."+s[1]+s[3:13]+'E'+"{:+03}".format(int(s[14:])+1)

        with zopen(file_name, "wt") as f:
            p = Poscar(self.structure)

            # use original name if it's been set (e.g. from Chgcar)
            comment = getattr(self, 'name', p.comment)

            lines = comment + "\n"
            lines += "   1.00000000000000\n"
            latt = self.structure.lattice.matrix
            lines += " %12.6f%12.6f%12.6f\n" % tuple(latt[0, :])
            lines += " %12.6f%12.6f%12.6f\n" % tuple(latt[1, :])
            lines += " %12.6f%12.6f%12.6f\n" % tuple(latt[2, :])
            if not vasp4_compatible:
                lines += "".join(["%5s" % s for s in p.site_symbols]) + "\n"
            lines += "".join(["%6d" % x for x in p.natoms]) + "\n"
            lines += "Direct\n"
            for site in self.structure:
                lines += "%10.6f%10.6f%10.6f\n" % tuple(site.frac_coords)
            lines += " \n"
            f.write(lines)
            a = self.dim

            def write_spin(data_type):
                lines = []
                count = 0
                f.write("   {}   {}   {}\n".format(a[0], a[1], a[2]))
                for (k, j, i) in itertools.product(list(range(a[2])),
                                                   list(range(a[1])),
                                                   list(range(a[0]))):
                    lines.append(_print_fortran_float(self.data[data_type][i, j, k]))
                    count += 1
                    if count % 5 == 0:
                        f.write(" " + "".join(lines) + "\n")
                        lines = []
                    else:
                        lines.append(" ")
                f.write(" " + "".join(lines) + " \n")
#                f.write("".join(self.data_aug.get(data_type, [])))
            write_spin(data_type)

def spin_density():
   proc_str="Reading Spin Charge Density From CHG File ..."
   procs(proc_str,1,sp='-->>')
   chg = CHGCAR.from_file("CHG")
   assert(chg.is_spin_polarized)
   proc_str="Writing Total Charge Density To Total.vasp File ..."
   procs(proc_str,2,sp='-->>')
   chg.write_file('Total.vasp','total')
   proc_str="Writing Spin Charge Density To Spin.vasp File ..."
   procs(proc_str,3,sp='-->>')
   chg.write_file('Spin.vasp','diff')

def charge_density_diff():
   print("input the name of charge density file like this")
   print("file_A file_B, or")
   print("it means rho(A)-rho(B))")
   print("file_A file_B file_C")
   print("it means rho(A)-rho(B)-rho(C)")
   wait_sep()
   in_str=""
   while in_str=="":
         in_str=input().strip()
   filenames=in_str.split()
   len_file=len(filenames)
   chg=[None]*len_file
   if len_file>=2:
      for i in range(len_file):
          fchg=filenames[i]
          proc_str="Reading Charge Density From "+fchg+" File ..."
          procs(proc_str,i+1,sp='-->>')
          chg[i] = CHGCAR.from_file(fchg)
      diff_chg=copy(chg[0])
      for j in range(len_file-1):
          diff_chg-=chg[j+1]
      proc_str="Writing Charge Density Difference to Diff.vasp  File ..."
      procs(proc_str,i+2,sp='-->>')
      diff_chg.write_file('Diff.vasp','total')
   else:
      print("you must input at least 2 files") 
      #charge_density_diff()
   
def spin_density_component():
   proc_str="Reading Spin Up/Down Density From CHG File ..."
   procs(proc_str,1,sp='-->>')
   chg = CHGCAR.from_file("CHG")
   assert(chg.is_spin_polarized)
   data=chg.spin_data.copy()
   chg.data={'spinup':data[Spin.up],'spindown':data[Spin.down]}
   proc_str="Writing Spin Up Density To Spinup.vasp File ..."
   procs(proc_str,2,sp='-->>')
   chg.write_file('Spinup.vasp','spinup')
   proc_str="Writing Spin Down Density To Spindown.vasp File ..."
   procs(proc_str,3,sp='-->>')
   chg.write_file('Spindown.vasp','spindown')

def chg_locp_average():
   chg=False
   print("which file would like to average: CHG or LOCPOT ?")
   wait_sep()
   in_str=""
   while in_str=="":
         in_str=input().strip()
   if in_str.lower()=='chg':
      filename='CHG'
      check_file(filename)
      grid_data = CHGCAR.from_file(filename)
      head_line="#%(key1)+s %(key2)+s"%{'key1':'Distance/Ang','key2':'Average Charge/(e/Ang**3)'}
      chg=True
   elif in_str.lower()=='locpot':
      filename='LOCPOT'
      check_file(filename)
      grid_data = Locpot.from_file(filename)
      head_line="#%(key1)+s %(key2)+s"%{'key1':'Distance/Ang','key2':'Average Potential/(eV)'}
   else:
      print('unknown file file: '+in_str)
      os._exit()

   step_count=1
   check_file(filename)
   proc_str="Reading Data From "+ filename +" File ..."
   procs(proc_str,step_count,sp='-->>')
   volume=grid_data.structure.volume 
   print("which direction would like to average: x y or z ?")
   wait_sep()
   in_str=""
   while in_str=="":
      in_str=input().strip().lower()
   if in_str=="x":
      selected_dir=0
   elif in_str=='y':
      selected_dir=1
   elif in_str=='z':
      selected_dir=2
   else:
       print("Unknow Direction!")
       return

   axis_grid=grid_data.get_axis_grid(selected_dir)
   if chg:
      aver_grid=grid_data.get_average_along_axis(selected_dir)/volume
   else:
      aver_grid=grid_data.get_average_along_axis(selected_dir)

   data=np.vstack((axis_grid,aver_grid)).T

   step_count+=1
   if chg:
      filename="average_CHG_"+in_str+'.dat'
   else:
      filename="average_LOCPOT_"+in_str+'.dat'

   proc_str="Writting Average Data to "+ filename +" File ..."
   procs(proc_str,step_count,sp='-->>')
   write_col_data(filename,data,head_line)
   
def optics_analysis():

   filename='vasprun.xml'
   step_count=1
   check_file(filename)
   proc_str="Reading Data From "+ filename +" File ..."
   procs(proc_str,step_count,sp='-->>')
   vsr=Vasprun(filename,parse_eigen=False)
   try:
     energy=np.array(vsr.dielectric[0])
     freq=energy/H
     real=np.array(vsr.dielectric[1])
     imag=np.array(vsr.dielectric[2])
   except:
     print('extracing data failed ')
     return

   head_line="#%(key1)+12s%(key2)+12s%(key3)+12s%(key4)+12s%(key5)+12s%(key6)+12s%(key7)+12s"%\
    {'key1':'Energy','key2':'xx','key3':'yy','key4':'zz','key5':'xy','key6':'yz','key7':'zx'}
   step_count+=1
   filename="AbsorbSpectrum.dat"
   proc_str="Writing Data To "+ filename +" File ..."
   procs(proc_str,step_count,sp='-->>')
   Absorb=np.zeros((len(freq),6))
   for i in range(6): 
       Absorb[:,i]=1.4142*freq*(np.sqrt(-real[:,i]+np.sqrt(imag[:,i]*imag[:,i]+real[:,i]*real[:,i])))/(C0*100)
   data=np.vstack((energy,Absorb.T)).T
   write_col_data(filename,data,head_line,e_fmt=True)

   step_count+=1
   filename="RefractiveSpectrum.dat"
   proc_str="Writing Data To "+ filename +" File ..."
   procs(proc_str,step_count,sp='-->>')
   Refrac=np.zeros((len(freq),6))
   for i in range(6): 
       Refrac[:,i]=(np.sqrt(real[:,i]+np.sqrt(imag[:,i]*imag[:,i]+real[:,i]*real[:,i])))/1.4142
   data=np.vstack((energy,Refrac.T)).T
   write_col_data(filename,data,head_line,e_fmt=True)

   step_count+=1
   filename="EnergyLossSpectrum.dat"
   proc_str="Writing Data To "+ filename +" File ..."
   procs(proc_str,step_count,sp='-->>')
   EnergyLoss=np.zeros((len(freq),6))
   for i in range(6): 
       EnergyLoss[:,i]=imag[:,i]/(imag[:,i]*imag[:,i]+real[:,i]*real[:,i])
   data=np.vstack((energy,EnergyLoss.T)).T
   write_col_data(filename,data,head_line,e_fmt=True)

   step_count+=1
   filename="ExtictionSpectrum.dat"
   proc_str="Writing Data To "+ filename +" File ..."
   procs(proc_str,step_count,sp='-->>')
   Extic=np.zeros((len(freq),6))
   for i in range(6): 
       Extic[:,i]=(np.sqrt(-real[:,i]+np.sqrt(imag[:,i]*imag[:,i]+real[:,i]*real[:,i])))/1.4142
   data=np.vstack((energy,Extic.T)).T
   write_col_data(filename,data,head_line,e_fmt=True)

   step_count+=1
   filename="ReflectivitySpectrum.dat"
   proc_str="Writing Data To "+ filename +" File ..."
   procs(proc_str,step_count,sp='-->>')
   Reflect=np.zeros((len(freq),6))
   for i in range(6): 
       Reflect[:,i]=((Refrac[:,i]-1)*(Refrac[:,i]-1)+Extic[:,i]*Extic[:,i])/((Refrac[:,i]+1)*(Refrac[:,i]+1)+Extic[:,i]*Extic[:,i])
   data=np.vstack((energy,Reflect.T)).T
   write_col_data(filename,data,head_line,e_fmt=True)

def aimd_analysis():
   pass

def elastic_analysis():

   filename='OUTCAR'
   step_count=1
   check_file(filename)
   proc_str="Reading Data From "+ filename +" File ..."
   procs(proc_str,step_count,sp='-->>')
   outcar=Outcar(filename)
   outcar.read_elastic_tensor()
   cij_tensor=np.array(outcar.data['elastic_tensor'])
   
   filename='vasprun.xml'
   step_count+=1
   check_file(filename)
   proc_str="Reading Data From "+ filename +" File ..."
   procs(proc_str,step_count,sp='-->>')
   vsr=Vasprun(filename)
   struct0=vsr.structures[0]
   natom=struct0.num_sites
   weight=struct0.composition.weight
   volume=struct0.lattice.volume

   ## converting the units
   volume *= 1.0e-30     ## from Angstrom to meters
   weight *= weight*1.0e-3 ## from gram to kg
   density = weight/(volume*Avogadro)

   asa=analyzer.SpacegroupAnalyzer(struct0)
#   lat_type=asa.get_crystal_system()
   
   crys_type=asa.get_lattice_type()

   ## Redefining the Cij matrix into the correct Voigt notation since VASP's OUTCAR has a different order
   ## In VASP: Columns and rows are listed as: 1, 2, 3, 6, 4, 5
   ## In this format OUTCAR's C44 values would be actually C66, C55 would be C44, and C66 would be C55. 
   ## OUTCAR follows the below order: 
   ## [C11 C12 C13 C16 C14 C15]
   ## [C21 C22 C23 C26 C24 C25]
   ## [C31 C32 C33 C36 C34 C35]
   ## [C61 C62 C63 C66 C64 C65]
   ## [C41 C42 C43 C46 C44 C45]
   ## [C51 C52 C53 C56 C54 C55] 

   cnew = np.zeros((6,6))
   snew = np.zeros((6,6))
   cnew = np.copy(cij_tensor)
   
   for j in range(0,6):
     cnew[3][j] = cij_tensor[4][j]
     cnew[4][j] = cij_tensor[5][j]
     cnew[5][j] = cij_tensor[3][j]
   
   ctemp=np.zeros((6,6))
   ctemp = np.copy(cnew)
   
   for i in range(0,6):
     cnew[i][3] = cnew[i][4]
     cnew[i][4] = cnew[i][5]
     cnew[i][5] = ctemp[i][3]
   
   # Change the units of Cij from kBar to GPa
   cnew=cnew/10.0;  
   proc_str="\n Modified elastic tensor in correct order (in GPa units)"
   print(proc_str)
   fmt="%7.3f "*6
   for i in range(6):
       print(fmt% tuple(cnew[i,:]))
    
   def check_symmetric(a, tol=1e-8):
       return np.allclose(a, a.T, atol=tol)

   print('\n Checking if the elastic tensor is symmetric: i.e. Cij = Cji:  %5s' % check_symmetric(cnew))
   print("\n Eigen Values of the elastic tensor:")
   evals = LA.eigvals(cnew)
   print(fmt% tuple(evals))
   if np.all(evals) > 0.0:
      print( "\n All eigen values are positive indicating elastic stability.")
   else:
      print("\n ATTENTION: One or more eigen values are negative indicating elastic instability.")
   calc_elastic_prop(cnew,snew,crys_type,density,weight,natom)


