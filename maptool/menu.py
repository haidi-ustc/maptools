#!/usr/bin/env python
from six.moves import input
from maptool.utils import sepline,box_center,wait_sep
from maptool.submenu import select_function
from maptool.version import logo,version
from maptool import NAME

def hello():
    box_center(ch='-',fill='-',sp="+")
    logo()
    box_center(ch='')
    version()
    box_center(ch='')
    box_center(ch='Written by Wang haidi')
    box_center(ch='URL  https://github.com/haidi-ustc')
    box_center(ch='Bug reports:(haidi@mail.ustc.edu.cn)')
    box_center(ch='-',fill='-',sp="+")

def byebye():
    box_center(ch='-',fill='-',sp="+")
    box_center(ch='*BYEBYE*')
    box_center(ch="Thanks for using "+NAME)
    box_center(ch="Have a nice day!!!")
    box_center(ch="Please cite: Nanoscale 9 (2), 850-855")
    box_center(ch="https://scholar.google.com/citations?hl=zh-CN&user=9PPScBEAAAAJ")
    box_center(ch='-',fill='-',sp="+")

def print_menu():
    sepline(ch=" structural operation ",sp='=')
    print('{} >>> {}'.format('a1','random operation'))
    print('{} >>> {}'.format('a2','covert operation'))
    print('{} >>> {}'.format('a3','build operation'))
    print('{} >>> {}'.format('a4','cleave operation'))
    print('{} >>> {}'.format('a5','strain operation'))
    print('{} >>> {}'.format('a6','2D structure operation'))

    sepline(ch=" structural analysis ",sp='=')
    print('{} >>> {}'.format('b1','structure symmetry'))
    print('{} >>> {}'.format('b2','structure finger print'))
    print('{} >>> {}'.format('b3','structure difference'))
    print('{} >>> {}'.format('b4','get primitive cell'))
    print('{} >>> {}'.format('b5','get conventional cell'))
    print('{} >>> {}'.format('b6','get XRD pattern'))

    sepline(ch=" vasp in/out tools ",sp='=')
    print('{} >>> {}'.format('c1','prepare input files'))
    print('{} >>> {}'.format('c2','analysis output files'))
#    print('{} >>> {}'.format('c3','summary output files'))

    #sepline(ch=" vasp auto calculation ",sp='=')
    #print('{} >>> {}'.format('d1','optimize structure'))
    #print('{} >>> {}'.format('d2','calculate band structure'))
    #print('{} >>> {}'.format('d3','calculate band structure HSE06'))
    #print('{} >>> {}'.format('d4','calculate dos'))
    #print('{} >>> {}'.format('d5','calculate dos by HSE06'))
    #print('{} >>> {}'.format('d6','calculate elastic properties'))
    #print('{} >>> {}'.format('d7','calculate phonon'))
    #print('{} >>> {}'.format('d8','execute MD simulation'))

    sepline(ch=" online retrieve ",sp='=')
    print('{} >>> {}'.format('e1','get band/dos by mp-ID'))
    print('{} >>> {}'.format('e2','get structure from materialsproject database'))
    print('{} >>> {}'.format('e3','get properties by mp-ID'))
    print('{} >>> {}'.format('e4','get phase graph'))
    sepline(ch="=",sp='=')

#TODO
#    sepline(ch=" siesta tools ",sp='=')
#    print('{} >>> {}'.format('si1','prepare input files'))
#    sepline(ch=" gaussian tools ",sp='=')
#    print('{} >>> {}'.format('gi1','prepare input files'))
#    sepline(ch=" nwchem tools ",sp='=')
#    print('{} >>> {}'.format('ni1','prepare input files'))
#    sepline(ch=" quantum espresso tools ",sp='=')
#    print('{} >>> {}'.format('qi1','prepare input files'))
#    sepline(ch=" lammps tools ",sp='=')
#    print('{} >>> {}'.format('li1','prepare input files'))

    print("")
    print("")
    print("88 >>> exit")
    print("your choice ?")
    wait_sep()
    in_str=""
    while in_str=="":
         in_str=input().strip()
    select_function(in_str)


if __name__=='__main__':
    hello()
    print_menu()
    byebye()

