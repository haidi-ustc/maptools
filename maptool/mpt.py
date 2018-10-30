#!/usr/bin/env python
from   __future__ import unicode_literals, print_function
import os
import sys
import string
import copy
from   datetime import datetime
from   optparse import OptionParser, OptionGroup
from   glob import glob
from   maptool.version import version,logo
from   maptool.menu import hello, print_menu,byebye
from   maptool.utils import box_center
from   maptool import info
from   time import time


__author__    = "Haidi Wang"
__copyright__ = "Copyright 2018"
__version__   = "0.1.1"
__maintainer__= "Haidi Wang"
__email__     = "haidi@mail.ustc.edu.cn"
__status__    = "Development"
__date__      = "May 16, 2018"


def main():
    parser = OptionParser()
    parser.set_defaults( Is_run = True,
                         struct_oper= False,
                         struct_ana= False,
                         vasp_inout=False,
                         vasp_auto=False,
                         online_extract=False,
                         cell_filename=None,
                         molecule_filename=None

                         )
    parser.add_option("--struct_operation",dest="struct_oper", action="store_true",help='structure operation')
    parser.add_option("--struct_analysis", dest="struct_ana", action="store_true",help='structure analysis')
    parser.add_option("--vasp_inout", dest="vasp_inout", action="store_true",help='create vasp input and analysis vasp output')
    parser.add_option("--vasp_auto", dest="vasp_auto", action="store_true",help='vasp automatic calculation')
    parser.add_option("-l","--logo", dest="Is_logo", action="store_true",help='show logo of pymat4vasp')
    parser.add_option("-v", "--version", dest="Is_version", action="store_true")
    parser.add_option("--clean", dest="Is_clean", action="store_true")
    parser.add_option("--ck", dest="check_file", action="store_true")
    parser.add_option("-r", "--run", dest="Is_run", action="store_true")
    parser.add_option("-i", "--info", dest="Is_info", action="store_true")
    parser.add_option("-c", "--cell", dest="cell_filename",action="store", type="string",help="Read crystal structure", metavar="FILE")
    parser.add_option("-m", "--molecule", dest="molecule_filename",action="store", type="string",help="Read molecule structure", metavar="FILE")

    (options, args) = parser.parse_args()

    rmvasp='''
rm -f IBZKPT CHG* CONTCAR DOSCAR EIGENVAL OSZICAR OUTCAR PCDAT XDATCAR WAVECAR  *err* vasprun.xml *out* *log* PARCHG* PROCAR ACF.dat AECCAR0 AECCAR1 AECCAR2 AVF.dat  BCF.dat ELFCAR BAND* REPORT DYNMAT FORCE_CONSTANTS *.conf COMPLETE RELAXED *.eps *.xyz vasp.std* *.json *.png
''' 
#    print options.Is_v,options.Is_clean,options.Is_plot,options.Is_job,options.Is_run
    ck_cmd="head POSCAR ; grep --color=auto VRHFIN POTCAR ; head KPOINTS"
    if options.check_file:
       options.Is_run=False
       os.system(ck_cmd)

    if options.Is_version:
        options.Is_run=False
        box_center(ch='-',fill='-',sp="+")
        version()
        box_center(ch='-',fill='-',sp="+")

    if options.Is_logo:
        options.Is_run=False
        box_center(ch='-',fill='-',sp="+")
        logo()
        box_center(ch='')
        box_center(ch='-',fill='-',sp="+")

    if options.Is_info:
       options.Is_run=False
       info()

    if options.Is_clean:
       options.Is_run=False
       os.system(rmvasp)
         
#    if options.Is_plot:
#       try:
#         import matplotlib.pyplot as pl
#         print("matplotlib is ok")
#       except:
#         print("Unable to load pyplot! Function PyGec.visual(...) will not work!")  

    if options.Is_run:
       T1=time()
       hello()
       print_menu()
       byebye()
       T2=time()
       print("Total Time: %.3f (s) "%(T2-T1))


if __name__=='__main__':
   main()

