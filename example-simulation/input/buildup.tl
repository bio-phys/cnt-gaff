source leaprc.lipid14
source leaprc.ff12SB
source leaprc.gaff
loadoff cnt.lib 
loadamberparams cnt.frcmod
loadamberparams frcmod.ionsjc_tip3p
system = loadpdb "system.pdb"
set system box {     83.224      83.456        84.1}
addions system Na+ 0
addions system Cl- 0
savepdb system system.pdb
saveamberparm system system.prmtop system.inpcrd
quit
