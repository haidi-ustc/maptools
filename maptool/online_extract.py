#!/usr/bin/env python
import os
from six.moves import input
from maptool.data_io import json_store
from maptool.vaspout import check_matplotlib
from maptool.utils import procs,sepline,wait_sep

from pymatgen import MPRester
from pymatgen.electronic_structure.plotter import DosPlotter,BSPlotter
from pymatgen.entries.compatibility import MaterialsProjectCompatibility
from pymatgen.entries.computed_entries import ComputedEntry
from pymatgen.analysis.phase_diagram import *


web="materials.org"

def check_apikey():
    try:
       apikey=os.environ['MAPI_KEY']
       mpr = MPRester(apikey)
       return(mpr)
    except:
       print(" you have to get a MAPI_KEY from "+web)
       print(" and execute following command:")
       print(" echo export MAPI_KEY='yourkey' >> ~/.bashrc")
       print(" source ~/.bashrc")
       os._exit()

def online_get_banddos():
    check_matplotlib()
    mpr=check_apikey()
    print("input the mp-ID")
    wait_sep()
    in_str=""
    while in_str=="":
          in_str=input().strip()
    mp_id=in_str
    step_count=1
    proc_str="Reading Data From "+ web +" ..."
    procs(proc_str,step_count,sp='-->>')
    data=mpr.get_entry_by_material_id(mp_id)
    sepline()
    print(data) 
    sepline()
    step_count+=1
    proc_str="Reading Band Data From "+ web +" ..."
    procs(proc_str,step_count,sp='-->>')
    band=mpr.get_bandstructure_by_material_id(mp_id)

    step_count+=1
    filename=mp_id+'_band.png'
    proc_str="Writing Data to "+filename +" File..."
    bsp=BSPlotter(band)
    procs(proc_str,step_count,sp='-->>')
    bsp.save_plot(filename=filename, img_format="png")
     
    step_count+=1
    proc_str="Reading DOS Data From "+ web +" ..."
    procs(proc_str,step_count,sp='-->>')
    dos=mpr.get_dos_by_material_id(mp_id)

    step_count+=1
    filename=mp_id+'_dos.png'
    proc_str="Writing Data to "+filename +" File..."
    dsp=DosPlotter()
    dsp.add_dos('Total',dos)
    procs(proc_str,step_count,sp='-->>')
    dsp.save_plot(filename=filename, img_format="png")

def online_get_structure():
    mpr=check_apikey()
    print('your choice ?')
    print('{} >>> {}'.format('1','get a structure by mp-ID'))
    print('{} >>> {}'.format('2','get a structure by fomular'))
    print('{} >>> {}'.format('3','get a structure by elements'))
    wait_sep()
    in_str=""
    while in_str=="":
          in_str=input().strip()
    choice=int(in_str)
    if choice==1:
       print("input the mp-ID")
       wait_sep()
       in_str=""
       while in_str=="":
             in_str=input().strip()
       mp_id=in_str
       struct=mpr.get_structure_by_material_id(mp_id)
       web="materials.org"
       step_count=1
       proc_str="Reading Data From "+ web +" ..."
       procs(proc_str,step_count,sp='-->>')
       filename=mp_id+'.vasp'
       step_count+=1
       proc_str="Writing Data to "+filename +" File..."
       procs(proc_str,step_count,sp='-->>')
       struct.to(filename=filename,fmt='POSCAR')
    elif choice==2:
       print("input the formula of structure")
       wait_sep()
       in_str=""
       while in_str=="":
             in_str=input().strip()
       formula=in_str
       mpid=mpr.query(criteria={'pretty_formula':formula},properties=['material_id'])
       web="materials.org"
       proc_str="Reading Data From "+ web +" ..."
       step_count=1
       procs(proc_str,step_count,sp='-->>')
       for i in range(len(mpid)):
            sid=mpid[i]['material_id']
            struct=mpr.get_structure_by_material_id(sid)
            filename=sid+'.cif'
            step_count+=1
            proc_str="Writing Data to "+filename +" File..."
            procs(proc_str,step_count,sp='-->>')
            struct.to(fmt='cif',filename=filename)
    elif choice==3:
       print("input the elements list")
       wait_sep()
       in_str=""
       while in_str=="":
             in_str=input().strip()
       elements=in_str.split()
       data=mpr.get_entries_in_chemsys(elements = elements)
       web="materials.org"
       proc_str="Reading Data From "+ web +" ..."
       step_count=1
       procs(proc_str,step_count,sp='-->>')
       for i in range(len(data)):
            if len(data[i].composition)==len(elements):
                sid=data[i].entry_id
                struct=mpr.get_structure_by_material_id(sid)
                filename=sid+'.cif'
                step_count+=1
                proc_str="Writing Data to "+filename +" File..."
                procs(proc_str,step_count,sp='-->>')
                struct.to(fmt='cif',filename=filename)
    else:
        print('unknow choice')

def online_get_phase_graph():
    mpr=check_apikey()
    compat=MaterialsProjectCompatibility()
    print("input the elements list")
    wait_sep()
    in_str=""
    while in_str=="":
          in_str=input().strip()
    elements=in_str.split()
    web="materials.org"
    proc_str="Reading Data From "+ web +" ..."
    step_count=1
    procs(proc_str,step_count,sp='-->>')
    unprocessed_entries=mpr.get_entries_in_chemsys(elements)
    processed_entries=compat.process_entries(unprocessed_entries)
    pd=PhaseDiagram(processed_entries)
    pdp=PDPlotter(pd,show_unstable=True)
    try:
       pdp.show()
    except:
       pass
    finally:
       step_count+=1
       filename='phase.png'
       proc_str="Writing Data to "+filename +" File..."
       procs(proc_str,step_count,sp='-->>')
       pdp.write_image(filename) 


def online_get_properties():
    mpr=check_apikey()
    print("input the mp-ID")
    wait_sep()
    in_str=""
    while in_str=="":
          in_str=input().strip()
    mp_id=in_str
    proc_str="Reading Data From "+ web +" ..."
    step_count=1
    procs(proc_str,step_count,sp='-->>')
    data=mpr.get_data(mp_id)
    filename=mp_id+'.json'
    proc_str="Writing Data To "+ filename +" File ..."
    step_count+=1
    procs(proc_str,step_count,sp='-->>')
    json_store(data,filename)

