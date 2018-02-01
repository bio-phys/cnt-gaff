#!/bin/bash

#TODO: Add documentation

for N in 08 ; do
	for LENGTH in 40 45 50; do
		for FUNC in none oh cooh coo; do
			# Make CNT structure
			./buildCstruct1_2.py -s armcnt -g $N $LENGTH --mol2 CNT_arm_$FUNC-$N-$LENGTH.mol2 -f $FUNC
			# Make CNT topology and Gromacs files
			./acpype.py -i CNT_arm_$FUNC-$N-$LENGTH.mol2 -c user > CNT_arm_$FUNC-$N-$LENGTH.log
		done
	done
done

exit 0
