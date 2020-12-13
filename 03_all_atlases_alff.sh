#!/bin/bash

sublist="/cbica/projects/alpraz_EI/scripts/alpraz_sublist.csv" 
measure="alffZ"
atlasdir=/cbica/projects/alpraz_EI/input/atlases/

for atlas in $(ls $atlasdir/*95.nii.gz); do
	fullname=$(basename $atlas)
	atlasName=${fullname%_thresh*}
	while IFS=, read subid sesid other; do
		for datavers in "GSR" "noGSR"; do
			(
			if [ $datavers = "GSR" ]; then
				outDir="/cbica/projects/alpraz_EI/data/SDK_TASK_GSR/${subid}/${sesid}"
				dataDir="/cbica/projects/alpraz_EI/data/SDK_TASK_GSR/xcpengine/sub-0${subid}/ses-00${sesid}/task-emotionid/space-MNI152NLin2009cAsym/alff/"
			else
				outDir="/cbica/projects/alpraz_EI/data/XCP_SDK_TASK_2020-03-25/${subid}/${sesid}"
				dataDir="/cbica/projects/alpraz_EI/data/XCP_SDK_TASK_2020-03-25/xcpengine/sub-0${subid}/ses-00${sesid}/task-emotionid/space-MNI152NLin2009cAsym/alff/"
			fi
			#Make the glasser timecourse
			if [ -f "${outDir}/${atlasName}_${measure}.1D" ]; then
				continue
			else
				3dROIstats -nomeanout -1DRformat -mask $atlas -nzmean ${dataDir}/sub-0${subid}_ses-00${sesid}_task-emotionid_space-MNI152NLin2009cAsym_${measure}.nii.gz > ${outDir}/${atlasName}_${measure}.1D
			fi
	
			)&
		done
		wait
	done <$sublist
done

