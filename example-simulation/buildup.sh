#!/bin/bash


module load gromacs/5.1.1-plumed2.2
module load amber/14 amberTools/14
module load vmd
module load python27/python/2.7


# Define the name of the carbon nanotube!

TUBENAME='CNT_arm_cooh-08-45'


# Get the carbon nanotube into its directory

mkdir nanotube
cp ../example-carbon-nanotubes/$TUBENAME.acpype/*_AC.* nanotube/
cp ../example-carbon-nanotubes/$TUBENAME.acpype/*_NEW.pdb nanotube/


# Create input directory

mkdir input
cd input


# Get all you need

cp ../nanotube/${TUBENAME}_AC.frcmod cnt.frcmod
cp ../nanotube/${TUBENAME}_AC.inpcrd cnt.inpcrd
cp ../nanotube/${TUBENAME}_AC.lib    cnt.lib
cp ../nanotube/${TUBENAME}_AC.prmtop cnt.prmtop
cp ../nanotube/${TUBENAME}_NEW.pdb   cnt-raw.pdb

cp ../charmm-gui/step5_assembly.pdb charmm-gui.pdb


# Convert charmm-gui structure to amber-readable structure

charmmlipid2amber.py -i charmm-gui.pdb  -c $AMBERHOME/AmberTools/src/etc/charmmlipid2amber/charmmlipid2amber.csv  -o bilayer.pdb


# Get the box size from the water coordinates

grep WAT bilayer.pdb | awk '{print $6, $7, $8}' > watercoords.txt
BOXAMB=$( python -c "import numpy as np; a = np.loadtxt('watercoords.txt'); d=0.1; print('%11.7g %11.7g %11.7g' % ( np.max(a[:,0])-np.min(a[:,0])+d, np.max(a[:,1])-np.min(a[:,1])+d, np.max(a[:,2
])-np.min(a[:,2])+d ) )" )
BOXGMX=$( python -c "import numpy as np; a = np.loadtxt('watercoords.txt'); d=0.1; print('%11.7g %11.7g %11.7g' % ( 0.1*(np.max(a[:,0])-np.min(a[:,0])+d), 0.1*(np.max(a[:,1])-np.min(a[:,1])+d), 0.1*(np.max(a[:,2
])-np.min(a[:,2])+d) ) )" )
rm watercoords.txt


# Put the carbon nanotube in the correct box

gmx editconf -f cnt-raw.pdb -o cnt-rot.pdb -rotate 90 0 0
gmx editconf -f cnt-rot.pdb -o cnt.pdb -box $BOXGMX -center 0 0 0


# Merge all together

grep 'ATOM\|TER' cnt.pdb      > all.pdb
echo 'TER'                   >> all.pdb
grep 'ATOM\|TER' bilayer.pdb >> all.pdb


# Remove molecules close to the nanotube)
# (Deletes TER, which causes errors in LEAP)

cat <<EOF > selection.vmd
mol new all.pdb
set sel [atomselect top "resname CNT or not same residue as (within 2.0 of resname CNT)"]
\$sel writepdb selected.pdb
exit
EOF
printf "source selection.vmd" | vmd
less selected.pdb | grep ATOM > out.pdb


# Correct for missing TER lines

cat <<EOF > addter.py
with open('out.pdb','r') as a:
    with open('system.pdb','w') as w:
        prev = '1'
        for line in a:
            l = line.split()
            if l[5] != prev:
                w.write("TER\n")
            w.write(line)
            prev = l[5]
EOF
python addter.py


# Create the input parameters using leap

cat << EOF > buildup.tl
source leaprc.lipid14
source leaprc.ff12SB
source leaprc.gaff
loadoff cnt.lib 
loadamberparams cnt.frcmod
loadamberparams frcmod.ionsjc_tip3p
system = loadpdb "system.pdb"
set system box {$BOXAMB}
addions system Na+ 0
addions system Cl- 0
savepdb system system.pdb
saveamberparm system system.prmtop system.inpcrd
quit
EOF
tleap -f buildup.tl
#less buildup.tl | xleap


# Leave directory
cd ..

exit
