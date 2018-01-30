#--------------------------
#  This project contains: 

 - buildCstruct 1.2 (Andrea Minoia, Martin Voegele) 
Generates atomistic structures of carbon nanotubes (CNTs) or graphite.
Can add three different functional groups to the rim (OH, COOH, COO-) of a CNT. 
Partial charges, if desired, are parameterized for use with the Amber force field. 
License: Free to use, modify and distribute

 - acpype.py (Alan Wilter Sousa da Silva)
Assigns generalized Amber (GAFF) parameters to organic molecules 
License: GNU General Public License V3.

 - two bash scripts (Martin Voegele)
Invoke both of the above to generate some example CNTs with various functional groups at the rim. 

 - a folder with ready-to-use example CNTs.


#-----------------
#  Requirements: 

 - Python 2.7 and Numpy

(for acpype.py:) 

 - Antechamber (from AmberTools preferably)

 - OpenBabel (optional, but strongly recommended)


#---------------
#  Literature:

 - M. Vögele, J. Köfinger, G. Hummer: 
   Simulations of Carbon Nanotube Porins in Lipid Bilayers.
   Faraday Discussions, 2018 (to be submitted)

(acpype)
 - A. W. Sousa da Silva, W. F. Vranken: 
   ACPYPE - AnteChamber PYthon Parser interfacE.
   BMC Research Notes 2012, 5:367 doi:10.1186/1756-0500-5-367
   http://www.biomedcentral.com/1756-0500/5/367

(Antechamber)
 - J. Wang, W. Wang, P. A. Kollman, D. A. Case: 
   Automatic atom type and bond type perception in molecular mechanical calculations. 
   Journal of Molecular Graphics and Modelling , 25, 2006, 247260.
 - J. Wang, R. M. Wolf, J. W. Caldwell, P. A. Kollman, D. A. Case:
   Development and testing of a general AMBER force field. 
   Journal of Computational Chemistry, 25, 2004, 1157-1174.


