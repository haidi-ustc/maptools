#!/usr/bin/env python

r'''
 This module used to control the convert the unit
'''
# Constant
Planck=6.6260755e-34    # Planck's constant, in Js
## Avogadro number
Avogadro = 6.022140857e+23 # Avogadro number
## Boltzmann constant
KB = 1.38064852e-23   # Boltzmann constant


C0 =2.99792458e+8
H  =4.13566733E-15

# Distance units
Bohr2Ang = 0.529177249  # Conversion of length from bohr to angstrom
Ang2Bohr = 1/Bohr2Ang
Ang2M= 1.0e-10 
M2Ang= 1.0/Ang2M 

# Energy units
Hartree2kCal = 627.5095 # Hartree to kcal/mol conversion
kCal2Hartree = 1/Hartree2kCal

eV2kCal = 23.061        # Conversion of energy in eV to energy in kcal/mol
kCal2eV = 1/eV2kCal

Hartree2Joule = 4.3597482e-18   # Hatree to Joule conversion factor
Joule2Hartree = 1/Hartree2Joule

eV2Hartree = Hartree2kCal/eV2kCal
Hartree2eV = 1/eV2Hartree

eV2Joule = 1.602e-19 
Joule2eV= 1/eV2Joule

Ry2eV =0.073498618 
eV2Py=1/Ry2eV

Ry2Hartree=0.50
Hartree2Ry=1/Ry2Hartree

#Pressure
Pa2GPa=1.0e-9

#others
eps_zero=1.0e-6
Len=70


