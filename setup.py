# coding: utf-8

import sys
import platform

from setuptools import setup, find_packages, Extension
from setuptools.command.build_ext import build_ext as _build_ext
default_prefix='maptool'

setup(
    name="maptool",
    packages=find_packages(),
    version="0.1.0",
    setup_requires=['pymatgen'],
    install_requires=["numpy>=1.9", "six", "requests", "ruamel.yaml>=0.15.6",
                      "monty>=0.9.6", "scipy>=1.0.0", "pydispatcher>=2.0.5",
                      "tabulate", "spglib>=1.9.9.44","pybtex>=0.21",'tqdm>=4.11.2',
                      "matplotlib>=1.5", "palettable>=2.1.1", "sympy", "pandas"],
#
#    extras_require={
#        ':python_version == "2.7"': [
#            'enum34',
#        ],
#        "provenance": ["pybtex"],
#        "ase": ["ase>=3.3"],
#        "vis": ["vtk>=6.0.0"],
#        "abinit": ["apscheduler==2.1.0"]},
    author="haidi",
    author_email="haidi@mail.ustc.edu.cn",
    maintainer="haidi",
    maintainer_email="haidi@mail.ustc.edu.cn",
    url="https://github.com/haidi-ustc/maptool",
    license="MIT",
    description=default_prefix+" is a postprocessing software for vasp, siesta, Gaussian, ABINIT and so on."
                "This software is based on Pymatgen and some open-source Python library.",
    keywords=["VASP", "gaussian", "ABINIT", "nwchem", "siesta","materials", "project","Quantum-Espresso"
              "electronic", "structure", "analysis", "phase", "diagrams"],
    classifiers=[
        "Programming Language :: Python :: 2",
        "Programming Language :: Python :: 2.7",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.5",
        "Programming Language :: Python :: 3.6",
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Topic :: Scientific/Engineering :: Information Analysis",
        "Topic :: Scientific/Engineering :: Physics",
        "Topic :: Scientific/Engineering :: Chemistry",
        "Topic :: Software Development :: Libraries :: Python Modules"
    ],
    entry_points={
          'console_scripts': [
              'mpt = maptool.mpt:main'],
 }
)
