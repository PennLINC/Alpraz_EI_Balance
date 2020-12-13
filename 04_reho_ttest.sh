#!/bin/bash
fd_thresh=.3
drug_list=""
placebo_list=""

datavers="acompcor" #"GSR"
measure="rehoZ"

/cbica/projects/alpraz_EI//data/XCP_SDK_TASK_2020-03-25/xcpengine/sub-014665/ses-002258/task-emotionid/space-MNI152NLin2009cAsym/alff/sub-014665_ses-002258_task-emotionid_space-MNI152NLin2009cAsym_alffZ.nii.gz 


while IFS=, read subid sesid day series group form drug other; do
	echo $subid $sesid
	if [ $datavers == "GSR" ]; then
		datadirname="XCP_SDK_TASK_GSR2020-03-30/"
	else
		datadirname="XCP_SDK_TASK_2020-03-25/"
	fi
	
	dataDir="/cbica/projects/alpraz_EI/data/${datadirname}/${subid}/${sesid}/"
	mkdir -p dataDir
	orgDir="/cbica/projects/alpraz_EI/data/${datadirname}/xcpengine/sub-0${subid}/ses-00${sesid}/task-emotionid/space-MNI152NLin2009cAsym/reho/"
	
	if [ $subid == "bblid" ]; then
		continue
	fi
	
	#Do a few checks to make sure data should be included
	#1. skip if we are missing either sessions
	if [ $(ls /cbica/projects/alpraz_EI//data/${datadir}/xcpengine/sub-0${subid}/ses-00*/task-emotionid/space-MNI152NLin2009cAsym/reho/sub*_task-emotionid_space-MNI152NLin2009cAsym_${measure}.nii.gz|wc -l) -lt 2 ];
	then 
		echo "$subid missing a file " 
		continue
	else
		if [ -f "$orgDir/sub-0${subid}_ses-00${sesid}_task-emotionid_space-MNI152NLin2009cAsym_${measure}.nii.gz" ]; then
			ln -sf "$orgDir/sub-0${subid}_ses-00${sesid}_task-emotionid_space-MNI152NLin2009cAsym_${measure}.nii.gz" ${dataDir}/${measure}.nii.gz
		fi
	fi
	#[[ ! -f ${dataDir}/amygdala_conn.nii.gz ]]&&continue
	
	#2. Check if there is too much motion in either scan
	fail=0
	for f in /cbica/projects/alpraz_EI/data/XCP_task-emotionid_acq-acompcor_2020-02-06/xcpengine/sub-0${subid}/ses-00*/task-emotionid/space-MNI152NLin2009cAsym/task/meanFD.txt; do
		fdtest=$(echo $(cat $f) ">" $fd_thresh|bc -l)
		[[ $fdtest = 1 ]]&&let "fail+=1" #mean FD is greater than the threshold
	done
	echo $fail
	if [ $fail -gt 0 ]; 
	then
		echo "$subid failed motion test"
		continue
	fi
	
	echo "$subid $sesid passed"

	# Now collect surviving scans.
	if [ $drug = 0 ]; then
		drug_list="$drug_list ${dataDir}/${measure}.nii.gz"
		echo "$subid $sesid added to drug list"
	elif [ $drug = 1 ]; then
		placebo_list="${placebo_list} ${dataDir}/${measure}.nii.gz"
		echo "$subid $sesid added to placebo list"
	fi
done </cbica/projects/alpraz_EI/scripts/alpraz_sub_groups.csv 

#export OMP_NUM_THREADS=40
3dttest++ -setA $drug_list -setB $placebo_list -paired -prefix /cbica/projects/alpraz_EI/output/${measure}_drug_ttest_${datavers}.nii.gz -overwrite -Clustsim 40 -mask /cbica/projects/alpraz_EI/input/coverage_mask_95perc.nii.gz