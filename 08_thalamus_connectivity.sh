#!/bin/bash
sublist="/cbica/projects/alpraz_EI/input/n1601_sublist.csv" 

while IFS=, read subid sesid other; do
	call="qsub -l short -N "thalamus_$subid_$sesid" proc_one_thalamus.sh $subid $sesid"
	eval $call

done <$sublist
