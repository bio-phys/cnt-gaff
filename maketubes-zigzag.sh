#!/bin/bash

#TODO: Add documentation

for N in 12; do
	for LENGTH in 40 45 50; do
		for FUNC in none oh cooh coo; do
			# Make CNT structure
			./buildCstruct1_2.py -s zigzagcnt -g $N $LENGTH --mol2 CNT_zigzag_$FUNC-$N-$LENGTH.mol2 -f $FUNC
			# Make CNT topology
			./acpype.py -i CNT_zigzag_$FUNC-$N-$LENGTH.mol2 -c user > CNT_zigzag_$FUNC-$N-$LENGTH.log
		done
	done
done

exit 0
