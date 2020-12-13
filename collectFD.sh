#!/bin/bash
outfile="/cbica/projects/alpraz_EI/input/alpraz_sublist_FD.csv"
fd_thresh=.5
echo "subid,sesid,day,series,group,form,drug,exists,FD,motion_pass">$outfile
while IFS=, read subid sesid day series group form drug other; do
	datadirname="XCP_task-emotionid_acq-acompcor_2020-02-06/"
	
	dataDir="/cbica/projects/alpraz_EI/data/${datadirname}/${subid}/${sesid}/"
	
	if [ $sesid = "fmriid" ]; then
		continue
	fi
	
	#Do a few checks to make sure data should be included
	#1. skip if we are missing either sessions
	if [ $(ls //cbica/projects/alpraz_EI/data/${datadirname}/${subid}/*/amygdala_conn_fisherz.nii.gz|wc -l) -lt 2 ];
	then 
		data_exists=0
		echo "$subid,$sesid,$day,$series,$group,$form,$drug,$data_exists,NaN,NaN">>$outfile
		continue
	else
		data_exists=1
	fi
	
	#2. Get FD
	fail=0
	for f in /cbica/projects/alpraz_EI/data/XCP_task-emotionid_acq-acompcor_2020-02-06/xcpengine/sub-0${subid}/ses-00*/task-emotionid/space-MNI152NLin2009cAsym/task/meanFD.txt; do
		fdtest=$(echo $(cat $f) ">" $fd_thresh|bc -l)
		[[ $fdtest = 1 ]]&&let "fail+=1" #mean FD is greater than the threshold
	done
	if [ $fail -gt 0 ]; 
	then
		motion_pass=0
	else
		motion_pass=1
	fi
	thisFD=$(cat /cbica/projects/alpraz_EI/data/XCP_task-emotionid_acq-acompcor_2020-02-06/xcpengine/sub-0${subid}/ses-00${sesid}/task-emotionid/space-MNI152NLin2009cAsym/task/meanFD.txt)
	
	echo "$subid,$sesid,$day,$series,$group,$form,$drug,$data_exists,$thisFD,$motion_pass">>$outfile

done </cbica/projects/alpraz_EI/scripts/alpraz_sub_groups.csv 
