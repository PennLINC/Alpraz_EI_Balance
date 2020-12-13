#!/bin/bash
fd_thresh=.3

datavers="acompcor" #"GSR"

echo "Subj drug emotion inputFile \\" > MVM_table.txt

while IFS=, read subid sesid day series group form drug other; do
	echo $subid $sesid
	if [ $datavers == "GSR" ]; then
		datadirname="XCP_SDK_TASK_GSR2020-03-30/"
	else
		datadirname="XCP_SDK_TASK_2020-03-25/"
	fi

	dataDir="/cbica/projects/alpraz_EI/data/${datadirname}/${subid}/${sesid}/"
	mkdir -p dataDir
	orgDir="/cbica/projects/alpraz_EI/data/${datadirname}/xcpengine/sub-0${subid}/ses-00${sesid}/task-emotionid/space-MNI152NLin2009cAsym/task/stats/"

	if [ $subid == "bblid" ]; then
		continue
	fi

	for measure in tstat1 tstat2 tstat3 tstat4 tstat5; do

		case $measure in
		tstat1)
			label="happy"
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
		if [ $(ls /cbica/projects/alpraz_EI/data/${datadirname}/xcpengine/sub-0${subid}/ses-00*/task-emotionid/space-MNI152NLin2009cAsym/task/stats/sub*_task-emotionid_space-MNI152NLin2009cAsym_${measure}.nii.gz|wc -l) -lt 2 ];
		then
			echo "$subid missing a file "
			continue
		else
			if [ -f "$orgDir/sub-0${subid}_ses-00${sesid}_task-emotionid_space-MNI152NLin2009cAsym_${measure}.nii.gz" ]; then
				ln -sf "$orgDir/sub-0${subid}_ses-00${sesid}_task-emotionid_space-MNI152NLin2009cAsym_${measure}.nii.gz" ${dataDir}/${label}_tstat.nii.gz
			fi
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
			#echo "$subid drug ${label} ${dataDir}/${label}_tstat.nii.gz" #just writing out to the screen. not important to do.
			echo "$subid drug ${label} ${dataDir}/${label}_tstat.nii.gz \\" >> MVM_table.txt #write to the MVM table inserting "drug" for the drug factor
		elif [ $drug = 1 ]; then
			#echo "$subid placebo ${label} ${dataDir}/${label}_tstat.nii.gz"
			echo "$subid placebo ${label} ${dataDir}/${label}_tstat.nii.gz \\" >> MVM_table.txt #inserting "placebo" for the drug factor
		fi
	done
done </cbica/projects/alpraz_EI/scripts/alpraz_sub_groups.csv
cat MVM_table.txt

3dMVM -prefix ${label}_mvm_${datavers}.nii.gz -jobs 20 \
	-mask /cbica/projects/alpraz_EI/input/coverage_mask_95perc.nii.gz \
	-wsVars "drug+emotion" \
	-num_glt 4 \
	-gltLabel 1 drug_by_happy-neutral_interaction -gltCode 1 'drug : 1*drug -1*placebo emotion : 1*happy -1*neutral' \
	-gltLabel 2 drug_by_angry-neutral_interaction -gltCode 2 'drug : 1*drug -1*placebo emotion : 1*angry -1*neutral' \
	-gltLabel 3 drug_by_sad-neutral_interaction -gltCode 3 'drug : 1*drug -1*placebo emotion : 1*sad -1*neutral' \
	-gltLabel 4 drug_by_fear-neutral_interaction -gltCode 4 'drug : 1*drug -1*placebo emotion : 1*fear -1*neutral' \
	-num_glf 5 \
	-glfLabel 1 happy_drug -glfCode 1 'emotion : 1*happy drug : 1*drug -1*placebo' \
	-glfLabel 2 fear_drug -glfCode 2 'emotion : 1*fear drug : 1*drug -1*placebo' \
	-glfLabel 3 sad_drug -glfCode 3 'emotion : 1*sad drug : 1*drug -1*placebo' \
	-glfLabel 4 angry_drug -glfCode 4 'emotion : 1*angry drug : 1*drug -1*placebo' \
	-glfLabel 5 neutral_drug -glfCode 5 'emotion : 1*neutral drug : 1*drug -1*placebo' \
	-dataTable @MVM_table.txt