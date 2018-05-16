#  This project contains: 

 - **buildCstruct 1.2** (Andrea Minoia, Martin Voegele)
Generates atomistic structures of carbon nanotubes (CNTs) or graphite.
Can add three different functional groups to the rim (OH, COOH, COO-) of a CNT. 
Partial charges, if desired, are parameterized for use with the Amber force field. 
See also http://chembytes.wikidot.com/buildcstruct.
License: Free to use, modify and distribute

 - **acpype.py** (Alan Wilter Sousa da Silva)
Assigns generalized Amber (GAFF) parameters to organic molecules.
This is included only for convenience and was taken unchanged as provided. 
See also http://www.ccpn.ac.uk/v2-software/software/ACPYPE-folder. 
License: GNU General Public License V3.

 - **Two example bash scripts** (Martin Voegele) and a folder with ready-to-use example CNTs:
The scripts invoke buildCstruct and acpype to generate some example CNTs with various functional groups at the rim. 
License: Free to use, modify and distribute


#  Requirements: 

 - Python 2.7 and Numpy

(for acpype.py:) 

 - Antechamber (from AmberTools preferably)

 - OpenBabel (optional, but strongly recommended)


# Installation

If all requirements are fulfilled (see above), no installation is required. 


# Get Started:

To build your CNT, invoke buildCstruct

    ./buildCstruct1_2.py -s [armcnt|zigzagcnt] -g [index n] [length in Angstrom] --mol2 [filename].mol2 -f [none|oh|cooh|coo]

To make the CNT topology, invoke acpype (for large systems, this might need a lot of memory)

    ./acpype.py -i [filename].mol2 -c user 

Both steps are scripted for various CNT geometries in
 -  maketubes-armchair.sh 
 -  maketubes-zigzag.sh


#  Literature:

 - M. Vögele, J. Köfinger, G. Hummer: 
   Simulations of Carbon Nanotube Porins in Lipid Bilayers.
   Faraday Discuss., 2018, Accepted Manuscript, DOI: 10.1039/C8FD00011E  
   http://pubs.rsc.org/en/content/articlelanding/2018/fd/c8fd00011e

(buildCstruct)
 - A. Minoia, L. Chen, D. Beljonne, L. Lazzaroni:
   Molecular Modeling Study of the Structure andStability of Polymer/Carbon Nanotube Interfaces.
   Polymer 53 (2012) 5480-5490

(acpype)
 - A. W. Sousa da Silva, W. F. Vranken: 
   ACPYPE - AnteChamber PYthon Parser interfacE.
   BMC Research Notes 2012, 5:367 doi:10.1186/1756-0500-5-367
   http://www.biomedcentral.com/1756-0500/5/367

(Antechamber)
 - J. Wang, W. Wang, P. A. Kollman, D. A. Case: 
   Automatic atom type and bond type perception in molecular mechanical calculations. 
   Journal of Molecular Graphics and Modelling, 25, 2006, 247260.
 - J. Wang, R. M. Wolf, J. W. Caldwell, P. A. Kollman, D. A. Case:
   Development and testing of a general AMBER force field. 
   Journal of Computational Chemistry, 25, 2004, 1157-1174.


