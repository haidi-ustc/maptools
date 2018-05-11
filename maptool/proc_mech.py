#!/usr/bin/env python
import numpy as np
from   maptool.constants import Avogadro,KB,Planck,Len

''' This script reads elastic tensor from OUTCAR, and calculates elastic and mechanical properties... 
 Authors: Sobhit Singh and Aldo Romero
          West Virginia University, Morgantown, USA
 modified by haidi
'''

def calc_elastic_prop(cnew,snew,crys_type,density,weight,natoms):
	elastic_const(cnew,snew,density,weight,natoms)
	
	crystal = np.array(['cubic', 'hexagonal', 'tetragonal', 'rhombohedral-1', 'rhombohedral-2', 'orthorhombic', 'monoclinic'])
	print("\n Supported system:\n ['cubic', 'hexagonal', 'tetragonal', 'rhombohedral-1',\n 'rhombohedral-2', 'orthorhombic', 'monoclinic']")
	print("\n")
	try:
           stability_test(cnew, crys_type)
	except:
           print("\n Choose one of the following: 'cubic', 'hexagonal', 'tetragonal', 'rhombohedral-1', 'rhombohedral-2', 'orthorhombic', 'monoclinic' \n")

def elastic_const(cnew, snew,density,weight,natoms):

	##Bulk: Voigt
	KV = (cnew[0][0] + cnew[1][1] + cnew[2][2]) + 2*(cnew[0][1] + cnew[1][2] + cnew[2][0])
	KV = KV/9.0

        ## Shear: Voigt
	GV = (cnew[0][0] + cnew[1][1] + cnew[2][2]) - (cnew[0][1] + cnew[1][2] + cnew[2][0]) + 3*(cnew[3][3] + cnew[4][4] + cnew[5][5])
	GV = GV/15.0

	## Young's: Voigt
	EV = (9*KV*GV)/(3*KV + GV)

	## Poisson's ratio: Voigt
	Nu_V = (3*KV - EV)/(6*KV)

	## P-wave modulus, M: Voigt
	MV = KV + (4*GV/3.0)

	#####  Reuss Method
	mc_new = np.matrix(cnew)
	mci_new = mc_new.I

	for i in range(0,6):
	    for j in range(0,6):
	        snew[i][j]=mci_new[i,j]

	## bulk: Reuss
	KR = (snew[0][0] + snew[1][1] + snew[2][2]) + 2*(snew[0][1] + snew[1][2] + snew[2][0])
	KR = 1.0/KR

	## Shear: Reuss
	GR = 4*(snew[0][0] + snew[1][1] + snew[2][2]) -4*(snew[0][1] + snew[1][2] + snew[2][0]) + 3*(snew[3][3] + snew[4][4] + snew[5][5])
	GR = 15.0/GR

	## Young's: Reuss
	ER = (9*KR*GR)/(3*KR + GR)

	## Poisson's ratio: Reuss
	Nu_R = (3*KR - ER)/(6*KR)

	## P-wave modulus, M: Reuss
	MR = KR + (4*GR/3.0)


	### Voigt-Reuss-Hill Approximation: average of both methods
	## VRH

	Kvrh = (KV + KR)/2.0
	Gvrh = (GV + GR)/2.0
	Evrh = (EV + ER)/2.0
	Nu_vrh = (Nu_V + Nu_R)/2.0
	Mvrh = (MV + MR)/2.0
	KG_ratio_V = KV/GV
	KG_ratio_R = KR/GR
	KG_ratio_vrh = Kvrh/Gvrh

	## Elastic Anisotropy
	## Zener anisotropy for cubic crystals only
	Az = 2*cnew[3][3]/(cnew[0][0] - cnew[0][1])

	##Ranganathan and Ostoja-Starzewski method: Phys. Rev. Lett. 101, 055504 (2008).
	## for any crystalline symmetry: Universal anisotropy index
	AU = (KV/KR) + 5*(GV/GR) - 6.0

	#""'Note: AU is a relative measure of anisotropy with respect to a limiting value. For example, AU does not prove that a crystal having AU = 3 has double the anisotropy of another crystal with AU = 1.5. I""'

	## log-Euclidean anisotropy parameter by Christopher M. Kube, AIP Advances 6, 095209 (2016)
	#""'AL  CV , CR   is based on the distance between the averaged stiffnesses CV and
	#CR , which is more appropriate. Clearly, AL  CV , CR   is zero when the crystallite is isotropic.

	AL = np.sqrt(5)*2.303*np.log(1 + (AU/5))
	

	print("\n \n                         Voigt     Reuss    Average")
	print(" "+'-'*(Len-1))
	print(" Bulk modulus   (GPa)  %9.3f %9.3f %9.3f " % (KV, KR, Kvrh))
	print(" Shear modulus  (GPa)  %9.3f %9.3f %9.3f " % (GV, GR, Gvrh))
	print(" Young modulus  (GPa)  %9.3f %9.3f %9.3f " % (EV, ER, Evrh))
	print(" Poisson ratio         %9.3f %9.3f %9.3f " % (Nu_V, Nu_R, Nu_vrh))
	print(" P-wave modulus  (GPa) %9.3f %9.3f %9.3f " % (MV, MR, Mvrh))
	print(" Bulk/Shear ratio      %9.3f %9.3f %9.3f (%s) " % (KG_ratio_V, KG_ratio_R, KG_ratio_vrh,  ductile_test(KG_ratio_vrh) ))
	print(" "+'-'*(Len-1))

	print("  \n \t \t Elastic Anisotropy \n ")
	print(" Zener anisotropy (true for cubic crystals only) Az = %10.5f" %Az)
	print(" Universal anisotropy index (Ranganathan and Ostoja-Starzewski method\n PRL 101, 055504 (2008)) Au = %10.5f" %AU)
	print(" Log-Euclidean anisotropy parameter by Christopher M. Kube\n AIP Advances 6, 095209 (2016) AL = %10.5f " %AL)


	## Calculation of elastic wave velocities and Debye temperature using the values obtained from Voigt-Reuss-Hill Approximation
	G = Gvrh*1.0e+9   ## converting from GPa to Pascal units (kg/ms^2)
	K = Kvrh*1.0e+9

	## transverse velocity: Navier's equation
	vt = np.sqrt((G/density))

        ## longitudinal velocity: Navier's equation
	vl = np.sqrt(((3*K + 4*G)/(3.0*density)))

        ## average 
	vm=1.0/(np.cbrt((1./3.)*(2./(vt*vt*vt)+1./(vl*vl*vl))))


	## Debye temperature calculated using  Orson Anderson's proposal [Ref- J. Phys. Chem. Solids (1963) 24, 909-917]
	q = np.sum(natoms)	
	theta = (Planck/KB)*vm*np.cbrt((3*q*Avogadro*density)/(4*(np.pi)*weight))

	## melting temperature using empirical relation from Ref: Johnston I, Keeler G, Rollins R and Spicklemire S 1996 
	##							  Solid State Physics Simulations, The Consortium for Upper-Level Physics Software (New York: Wiley)
	Tm = 607 + 9.3*Kvrh

	print(" Elastic wave velocities calculated using Navier's equation  (in m/s units) \n")
	print(" "+'-'*(Len-1))
	print(" Longitudinal wave velocity (vl) : %10.5f " %vl)
	print(" Transverse wave velocity (vt) : %10.5f " %vt)
	print(" Average wave velocity (vm) : %10.5f " %vm)
	print(" Debye temperature  (in K) : %10.5f " %theta)
	print(" \n Melting temperature calculated from empirical relation:\n Tm = 607 + 9.3*Kvrh \pm 555 (in K)")
	print(" Tm (in K)=  %10.5f (plus-minus 555 K) " % Tm)
	print(" "+'-'*(Len-1))


	

def ductile_test(ratio):
	if(ratio > 1.75):
		return "ductile"
	else:
		return "brittle"


def stability_test(matrix, crystaltype):
 c = np.copy(matrix)

 if(crystaltype =="cubic"):
   print(" Cubic crystal system \n")
   print(" Born stability criteria for the stability of cubic system are :\n [Ref- Mouhat and Coudert, PRB 90, 224104 (2014)]  \n")
   print(" (i) C11 - C12 > 0;    (ii) C11 + 2C12 > 0;   (iii) C44 > 0 \n ")

   ## check (i)   keep in mind list starts with 0, so c11 is stored as c00
   if(c[0][0] - c[0][1] > 0.0):
      print(" Condition (i) satified.")
   else:
      print(" Condition (i) NOT satisfied.")

   if(c[0][0] + 2*c[0][1] > 0.0):
      print(" Condition (ii) satified.")
   else:
      print(" Condition (ii) NOT satisfied.")

   if(c[3][3] > 0.0):
      print(" Condition (iii) satified.")
   else:
      print(" Condition (iii) NOT satisfied.")


 if(crystaltype =="hexagonal"):
   print(" Hexagonal crystal system \n")
   print(" Born stability criteria for the stability of hexagonal system are :\n [Ref- Mouhat and Coudert, PRB 90, 224104 (2014)]  \n")
   print(" (i) C11 - C12 > 0;    (ii) 2*C13^2 < C33(C11 + C12);   (iii) C44 > 0 \n ")

   ## check (i)   keep in mind list starts with 0, so c11 is stored as c00
   if(c[0][0] - c[0][1] > 0.0):
      print(" Condition (i) satified.")
   else:
      print(" Condition (i) NOT satisfied.")

   if(2*(c[0][2]*c[0][2]) < c[2][2]*(c[0][0] + c[0][1])):
      print(" Condition (ii) satified.")
   else:
      print(" Condition (ii) NOT satisfied.")

   if(c[3][3] > 0.0):
      print(" Condition (iii) satified.")
   else:
      print(" Condition (iii) NOT satisfied.")


 if(crystaltype =="tetragonal"):
   print(" Tetragonal crystal system \n")
   print(" Born stability criteria for the stability of Tetragonal system are :\n [Ref- Mouhat and Coudert, PRB 90, 224104 (2014)]  \n")
   print(" (i) C11 - C12 > 0;    (ii) 2*C13^2 < C33(C11 + C12);   (iii) C44 > 0; \n (iv) C66 > 0;    (v) 2C16^2 < C66*(C11-C12) \n ")

   ## check (i)   keep in mind list starts with 0, so c11 is stored as c00
   if(c[0][0] - c[0][1] > 0.0):
      print(" Condition (i) is satified.")
   else:
      print(" Condition (i) is NOT satisfied.")

   if(2*(c[0][2]*c[0][2]) < c[2][2]*(c[0][0] + c[0][1])):
      print(" Condition (ii) is satified.")
   else:
      print(" Condition (ii) is NOT satisfied.")

   if(c[3][3] > 0.0):
      print(" Condition (iii) is satified.")
   else:
      print(" Condition (iii) is NOT satisfied.")

   if(c[5][5] > 0.0):
      print(" Condition (iv) is satified.")
   else:
      print(" Condition (iv) is NOT satisfied.")

   if(2*c[0][5]*c[0][5] < c[5][5]*(c[0][0] - c[0][1])):
      print(" Condition (v) is satified.")
   else:
      print(" Condition (v) is NOT satisfied.")


 if(crystaltype =="rhombohedral-1"):
   print(" Rhombohedral (class-1): i.e. structures with point group: 3m, -3m and 32 \n")
   print(" Born stability criteria for the stability of Rhombohedral-1 class system are :\n [Ref- Mouhat and Coudert, PRB 90, 224104 (2014)]  \n")
   print(" (i) C11 - C12 > 0;    (ii) C13^2 < (1/2)*C33(C11 + C12);\n (iii) C14^2 < (1/2)*C44*(C11-C12) = C44*C66;   (iv)  C44 > 0; \n ")

   ## check (i)   keep in mind list starts with 0, so c11 is stored as c00
   if(c[0][0] - c[0][1] > 0.0):
      print(" Condition (i) is satified.")
   else:
      print(" Condition (i) is NOT satisfied.")

   if((c[0][2]*c[0][2]) < (0.5)*c[2][2]*(c[0][0] + c[0][1])):
      print(" Condition (ii) is satified.")
   else:
      print(" Condition (ii) is NOT satisfied.")

   if(c[0][3]*c[0][3] < 0.5*c[3][3]*(c[0][0] - c[0][1])):
      print(" Condition (iii) is satified.")
   else:
      print(" Condition (iii) is NOT satisfied.")

   if(c[3][3] > 0.0):
      print(" Condition (iv) is satified.")
   else:
      print(" Condition (iv) is NOT satisfied.")


 if(crystaltype =="rhombohedral-2"):
   print(" Rhombohedral (class-2): i.e structures with point group: 3, -3 \n")
   print(" Born stability criteria for the stability of Rhombohedral-1 class system are :\n [Ref- Mouhat and Coudert, PRB 90, 224104 (2014)]  \n")
   print(" (i) C11 - C12 > 0;    (ii) C13^2 < (1/2)*C33(C11 + C12); \n (iii) C14^2 + C15^2 < (1/2)*C44*(C11-C12) = C44*C66;   (iv)  C44 > 0;  Note: C15 is added.. \n ")

   ## check (i)   keep in mind list starts with 0, so c11 is stored as c00
   if(c[0][0] - c[0][1] > 0.0):
      print(" Condition (i) is satified.")
   else:
      print(" Condition (i) is NOT satisfied.")

   if(c[0][2]*c[0][2] < (0.5)*c[2][2]*(c[0][0] + c[0][1])):
      print(" Condition (ii) is satified.")
   else:
      print(" Condition (ii) is NOT satisfied.")

   if(c[0][3]*c[0][3] + c[0][4]*c[0][4]  < 0.5*c[3][3]*(c[0][0] - c[0][1])):
      print(" Condition (iii) is satified.")
   else:
      print(" Condition (iii) is NOT satisfied.")

   if(c[3][3] > 0.0):
      print(" Condition (iv) is satified.")
   else:
      print(" Condition (iv) is NOT satisfied.")


 if(crystaltype =="orthorhombic"):
   print(" Orthorhombic crystal system.... \n")
   print(" Born stability criteria for the stability of Orthorhombic systems are :\n [Ref- Mouhat and Coudert, PRB 90, 224104 (2014)]  \n")
   print(" (i) C11 > 0;   (ii) C11*C22 > C12^2; \n (iii) C11*C22*C33 + 2C12*C13*C23 - C11*C23^2 - C22*C13^2 - C33*C12^2 > 0; \n (iv)  C44 > 0;   (v)  C55 > 0 ;   (vi)  C66 > 0 \n")
   ## check (i)   keep in mind list starts with 0, so c11 is stored as c00
   if(c[0][0]  > 0.0):
      print(" Condition (i) is satified.")
   else:
      print(" Condition (i) is NOT satisfied.")

   if(c[0][0]*c[1][1] > c[0][1]*(c[0][1])):
      print(" Condition (ii) is satified.")
   else:
      print(" Condition (ii) is NOT satisfied.")

   if(c[0][0]*c[1][1]*c[2][2] + 2*c[0][1]*c[0][2]*c[1][2] - c[0][0]*c[1][2]*c[1][2] - c[1][1]*c[0][2]*c[0][2] - c[2][2]*c[0][1]*c[0][1] > 0 ):
      print(" Condition (iii) is satified.")
   else:
      print(" Condition (iii) is NOT satisfied.")

   if(c[3][3] > 0.0):
      print(" Condition (iv) is satified.")
   else:
      print(" Condition (iv) is NOT satisfied.")

   if(c[4][4] > 0.0):
      print(" Condition (iv) is satified.")
   else:
      print(" Condition (iv) is NOT satisfied.")

   if(c[5][5] > 0.0):
      print(" Condition (iv) is satified.")
   else:
      print(" Condition (iv) is NOT satisfied.")


 if(crystaltype =="monoclinic"):
   print("Independent elastic contant: C11 C22 C33 C44 C55 C66 C12 C13 C23 C15 C25 C35 C46")
   print(c[0][0],c[1][1],c[2][2],c[3][3],c[4][4],c[5][5])
   print(c[0][1],c[0][2],c[1][2])
   print(c[0][4],c[1][4],c[2][4])
   print(c[3][5])
   print(" Monoclinic crystal system.... \n")
   print(" Born stability criteria for the stability of monoclinic systems are :\n [Ref- Mouhat and Coudert, PRB 90, 224104 (2014), and Wu et al. PRB 76, 054115 (2007)]  \n")
   print(" (i) C11 > 0;  (ii)  C22 > 0; (iii)  C33 > 0;\n (iv)  C44 > 0;   (v)  C55 > 0 ;   (vi)  C66 > 0  ")
   print(" (vii) [C11 + C22 + C33 + 2*(C12 + C13 + C23)] > 0;  \n (viii)  C33*C55 - C35^2 > 0;   (ix)  C44*C66 - C46^2 > 0;   (x) C22 + C33 - 2*C23  > 0 ")
   print(" (xi) C22*(C33*C55 - C35^2) + 2*C23*C25*C35 - (C23^2)*C55 - (C25^2)*C33   > 0  ")
   print(" (xii)  2*[C15*C25*(C33*C12 - C13*C23) + C15*C35*(C22*C13 - C12*C23) + C25*C35*(C11*C23 - C12*C13)] - [C15*C15*(C22*C33 - C23^2) + C25*C25*(C11*C33 - C13^2) + C35*C35*(C11*C22 - C12^2)] + C55*g > 0  ")
   print("          where, g = [C11*C22*C33 - C11*C23*C23 - C22*C13*C13 - C33*C12*C12 + 2*C12*C13*C23 ] ")

   for i in range(0, 6):
      if(c[i][i]  > 0.0):
        print(" Condition (%2d) is satified." % (i+1))
      else:
        print(" Condition (%2d) is NOT satified." % (i+1))


   if(c[0][0] + c[1][1] + c[2][2] + 2*(c[0][1] + c[0][2] + c[1][2]) > 0 ):
      print(" Condition (vii) is satified.")
   else:
      print(" Condition (vii) is NOT satisfied.")

   if (c[2][2]*c[4][4] - c[2][4]*c[2][4] > 0):
      print(" Condition (viii) is satified.")
   else:
      print(" Condition (viii) is NOT satisfied.")


   if(c[3][3]*c[5][5] - c[3][5]*c[3][5] > 0.0):
      print(" Condition (ix) is satified.")
   else:
      print(" Condition (ix) is NOT satisfied.")


   if(c[1][1] + c[2][2] - 2*c[1][2] > 0.0):
      print(" Condition (x) is satified.")
   else:
      print(" Condition (x) is NOT satisfied.")


   if(c[1][1]*(c[2][2]*c[2][4] - c[2][4]*c[2][4]) + 2*c[1][2]*c[1][4]*c[2][4] - c[1][4]*c[1][4]*c[2][2] > 0.0):
      print(" Condition (xi) is satified.")
   else:
      print(" Condition (xi) is NOT satisfied.")

   g = (c[0][0]*c[1][1]*c[2][2]) - (c[0][0]*c[1][2]*c[1][2]) - (c[1][1]*c[0][2]*c[0][2]) - (c[2][2]*c[0][1]*c[0][1]) + 2.0*(c[0][1]*c[0][2]*c[1][2])
   h1 =  2*(c[0][4]*c[1][4]*(c[2][2]*c[0][1] - c[0][2]*c[1][2]) + c[0][4]*c[2][4]*(c[1][1]*c[0][2] - c[0][1]*c[1][2]) + c[1][4]*c[2][4]*(c[0][0]*c[1][2] - c[0][1]*c[0][2])) 
   h2 = (c[0][4]*c[0][4]*(c[1][1]*c[2][2] - c[1][2]*c[1][2]) + c[1][4]*c[1][4]*(c[0][0]*c[2][2] - c[0][2]*c[0][2]) + c[2][4]*c[2][4]*(c[0][0]*c[1][1] - c[0][1]*c[0][1]))
   x = h1 - h2 + c[4][4]*g


   if(x > 0.0):
      print(" Condition (xii) is satified.")
   else:
      print(" Condition (xii) is NOT satisfied.")

