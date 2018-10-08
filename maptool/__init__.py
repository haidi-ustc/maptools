#!/usr/bin/env python
from   __future__ import unicode_literals, print_function
import logging
import os

__author__    = "Haidi Wang"
__copyright__ = "Copyright 2018"
__version__   = "0.1.0"
__maintainer__= "Haidi Wang"
__email__     = "haidi@mail.ustc.edu.cn"
__status__    = "Development"
__date__      = "May 16, 2018"


DEBUG = True
NAME  = 'maptool'


mpt_log = logging.getLogger(__name__)
mpt_log.setLevel(logging.DEBUG)
mpt_logf = logging.FileHandler(os.getcwd()+os.sep+'mpt.log')
mpt_logf_formatter=logging.Formatter('%(asctime)s - %(levelname)s : %(message)s')
#mpt_logf_formatter=logging.Formatter('%(asctime)s - %(name)s - [%(filename)s:%(funcName)s - %(lineno)d ] - %(levelname)s \n %(message)s')
mpt_logf.setFormatter(mpt_logf_formatter)
mpt_log.addHandler(mpt_logf)

try:
   import skued
   SKUED=True
except:
   SKUED=False

try: 
   import ase
   ASE=True
   print(ASE)
except:
   ASE=False

def info():
    """
        Show basic information about maptool, its location and version.
        Also information about other libraries used by maptool
        both mandatory and optional
    """

    print('Maptool\n--------\n')
    print('Version: ' + __version__)
    print('Path:    ' + __path__[0])
    print('Date:    ' + __date__)
    print()

    import sys
    print('Python version=' + sys.version + '\n')

    try:
        mm = __import__('pymongo')
        print('%10s %10s   %s' % ('pymongo', mm.version, mm.__path__[0]))
    except ImportError:
        print('pymongo Not Found')

    for modui in ['numpy', 'scipy', 'mayavi', 'matplotlib', 'tqdm',
                  'future', 'nose', 'coverage', 'spglib', 'pyhull', 'pymatgen', 'qmpy', ]:
        try:
            mm = __import__(modui)
            print('%10s %10s   %s' % (modui, mm.__version__, mm.__path__[0]))
        except ImportError:
            print('%10s %10s Not Found' % (modui, ''))

    if ASE:
        import ase
        #from ase import version as ase_version
        print('%10s %10s   %s' % ('ase', ase.__version__, ase.__path__[0]))
    else:
        print('%10s %10s Not Found' % ('ase', ''))
