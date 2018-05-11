#!/usr/bin/env python
from    __future__ import division, unicode_literals, print_function
import os
import json
import numpy as np
from   monty.io import zopen


def json_store(data,filename):
    with open(filename, 'w') as json_file:
        json_file.write(json.dumps(data,indent=4))

def json_read(filename):
    jsObj = json.load(open(filename)) 
    return jsObj

def write_col_data(file_name,data,head_line='',nspace=None,str_data=False,e_fmt=False,sp_fmt=None):

   # directly write the string data
   if str_data:
      with zopen(file_name, "wt") as f:
         f.write(data)
      return

   # set up the output format
   row,col=np.shape(data)
   if sp_fmt is None:
      if e_fmt:
         fmt="    %8.3e"*col+'\n'
      else:
         fmt="%12.6f"*col+'\n'
   else:
      fmt=sp_fmt

   # write data without linefeed
   if nspace is None:
      with zopen(file_name, "wt") as f:
          lines = head_line + "\n"
          for idata in data:
              lines += fmt % tuple(idata)
          lines += " \n"
          f.write(lines)

   # write data with linefeed
   elif isinstance(nspace,int):
      lcount=0
      with zopen(file_name, "wt") as f:
          lines = head_line + "\n"
          for idata in data:
              lines += fmt % tuple(idata)
              lcount+=1
              if lcount%nspace==0:
                 lines+='\n'
          lines += " \n"
          f.write(lines)
   else:
       print("unknow format ")
       os._exit()
