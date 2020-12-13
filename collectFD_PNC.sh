#!/bin/bash
outfile="/cbica/projects/alpraz_EI/input/PNC_sublist_FD.csv"
fd_thresh=.3
sublist="/cbica/projects/alpraz_EI/input/n1601_sublist.csv" 
echo "subid,sesid,exists,FD,motion_pass">$outfile
while IFS=, read subid sesid other; do
	
	dataDir="/cbica/projects/alpraz_EI/data/PNC/IDEMO_ACOMPCOR_GSR/${subid}/${sesid}/"
	
	if [ $sesid = "scanid" ]; then
		continue
	fi
	
	#Do a few checks to make sure data should be included
	#1. skip if we are missing data
	if [ ! -f ${dataDir}/schaefer200x7_aal_CC_000.netcc ];
	then 
		data_exists=0
		echo "$subid,$sesid,$data_exists,NaN,NaN">>$outfile
		continue
	else
		data_exists=1
	fi
	
	#2. Get FD
	measure=$(cat "/cbica/projects/alpraz_EI/data/PNC/IDEMO_ACOMPCOR_GSR/xcpengine/sub-${subid}/ses-PNC1/task-idemo/space-MNI152NLin2009cAsym/sub-${subid}_ses-PNC1_task-idemo_space-MNI152NLin2009cAsym_quality.csv" |cut -d "," -f8|sed 2d)
	FD=$(cat "/cbica/projects/alpraz_EI/data/PNC/IDEMO_ACOMPCOR_GSR/xcpengine/sub-${subid}/ses-PNC1/task-idemo/space-MNI152NLin2009cAsym/sub-${subid}_ses-PNC1_task-idemo_space-MNI152NLin2009cAsym_quality.csv" |cut -d "," -f8|sed 1d)
	fdtest=$(echo $FD "<" $fd_thresh|bc -l)

	echo "$subid,$sesid,$data_exists,$FD,$fdtest">>$outfile

done <$sublist
