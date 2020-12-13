#!/bin/bash
set -xe
datadir=/cbica/projects/alpraz_EI/input/atlases/
mkdir -p datadir

threshold=.95

for atlasPath in $(find /cbica/projects/alpraz_EI/data/XCP_SDK_TASK_GSR2020-03-30/xcpengine/sub-013372/ses-001539/task-emotionid/space-MNI152NLin2009cAsym/roiquant/* -maxdepth 1 -type d); do
	atlasName=${atlasPath##*/}
	atlas=${atlasPath}/*${atlasName}.nii.gz
	echo $atlasName
	echo $atlas
	
	cp $atlas $datadir/$atlasName.nii.gz

	# Get mean coverage for all parcels in the atlas
	3dROIstats -mask $atlas -nobriklab -1DRformat  ${datadir}/mask_mean.nii.gz > ${datadir}/atlas_means.txt

	# Run Rscript to threshold the atlas
	Rscript CreateThresholdedAtlas.R $atlas $threshold
	
	mv /cbica/projects/alpraz_EI/input/atlases/atlas_threshold_0${threshold}.nii.gz /cbica/projects/alpraz_EI/input/atlases/${atlasName}_threshold_0${threshold}.nii.gz
done