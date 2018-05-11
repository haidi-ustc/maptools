#!/usr/bin/env python
from maptool.constants import Len
from maptool.utils import box_center
import numpy as np 
r"""
logo and version description
"""

def logo():
     if np.random.rand()>0.5:
        s=[None]*4
        s[0]="  _   _   _   _   _   _   _  "
        s[1]=" / \ / \ / \ / \ / \ / \ / \ "
        s[2]="( m | a | p | t | o | o | l )"
        s[3]=" \_/ \_/ \_/ \_/ \_/ \_/ \_/ "
     else:
        s=[None]*6
        s[0]="                       _              _ "
        s[1]=" _ __ ___   __ _ _ __ | |_ ___   ___ | |"
        s[2]="| '_ ` _ \ / _` | '_ \| __/ _ \ / _ \| |"
        s[3]="| | | | | | (_| | |_) | || (_) | (_) | |"
        s[4]="|_| |_| |_|\__,_| .__/ \__\___/ \___/|_|"
        s[5]="                |_|                     "
              
     for i in range(len(s)):
         box_center(ch=s[i])
def version():
     version=" VERSION 0.1.0 (2018.2.26) "
     box_center(ch=version)
