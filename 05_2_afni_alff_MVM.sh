#!/bin/bash
fd_thresh=.3
drug_list=""
placebo_list=""

datavers="acompcor" #"GSR"

echo "Subj drug inputFile \\" > MVM_alff_table.txt

while IFS=, read subid sesid day series group form drug other; do
	echo $subid $sesid
	if [ $datavers == "GSR" ]; then
		datadirname="XCP_SDK_TASK_GSR2020-03-30/"
	else
		datadirname="XCP_SDK_TASK_2020-03-25/"
	fi

	dataDir="/cbica/projects/alpraz_EI/data/${datadirname}/${subid}/${sesid}/"

	if [ $subid == "bblid" ]; then
		continue
	fi

	for measure in alffZ; do

		case $measure in
		alffZ)
			label="alffZ"
			;;
		tstat2)
			label="fear"
			;;
		tstat3)
			label="angry"
			;;
		tstat4)
			label="neutral"
			;;
		tstat5)
			label="sad"
			;;
		esac

		#Do a few checks to make sure data should be included
		#1. skip if we are missing either sessions
		if [ $(ls /cbica/projects/alpraz_EI/data/${datadirname}/${subid}/*/${measure}.nii.gz|wc -l) -lt 2 ];
		then
			echo "$subid missing a file "
			continue
		fi
		#[[ ! -f ${dataDir}/amygdala_conn.nii.gz ]]&&continue

		#2. Check if there is too much motion in either scan
		fail=0
		for f in /cbica/projects/alpraz_EI/data/XCP_task-emotionid_acq-acompcor_2020-02-06/xcpengine/sub-0${subid}/ses-00*/task-emotionid/space-MNI152NLin2009cAsym/task/meanFD.txt; do
			fdtest=$(echo $(cat $f) ">" $fd_thresh|bc -l)
			[[ $fdtest = 1 ]]&&let "fail+=1" #mean FD is greater than the threshold
		done
		if [ $fail -gt 0 ];
		then
			continue
		fi


		# Now collect surviving scans.
		if [ $drug = 0 ]; then
			echo "$subid drug ${dataDir}/${measure}.nii.gz"
			echo "$subid drug ${dataDir}/${measure}.nii.gz \\" >> MVM_alff_table.txt
		elif [ $drug = 1 ]; then
			echo "$subid placebo ${dataDir}/${measure}.nii.gz"
			echo "$subid placebo ${dataDir}/${measure}.nii.gz \\" >> MVM_alff_table.txt
		fi
	done
done </cbica/projects/alpraz_EI/scripts/alpraz_sub_groups.csv
cat MVM_table.txt

3dMVM -prefix ${label}_mvm_${measure}_${datavers}.nii.gz -jobs 20 \
	-mask /cbica/projects/alpraz_EI/input/coverage_mask_95perc.nii.gz \
	-wsVars "drug" \
	-dataTable @MVM_alff_table.txt