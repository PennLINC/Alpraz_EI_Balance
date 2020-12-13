#!/bin/bash

sublist="/cbica/projects/alpraz_EI/scripts/alpraz_sublist.csv" 
#!/bin/bash

sublist="/cbica/projects/alpraz_EI/input/n1601_sublist.csv" 



while IFS=, read subid sesid other; do
	outDir="/cbica/projects/alpraz_EI/data/PNC/IDEMO_ACOMPCOR_GSR/${subid}/${sesid}"
	dataDir="/cbica/projects/alpraz_EI/data/PNC/IDEMO_ACOMPCOR_GSR/xcpengine/sub-${subid}/ses-PNC1/task-idemo/space-MNI152NLin2009cAsym/task/stats/"
	dataFile="${dataDir}/sub-${subid}_ses-PNC1_task-idemo_space-MNI152NLin2009cAsym_res4d.nii.gz"
	datavers="GSR"

	#Make the thalamus timecourses
	# 3dcalc -a ${dataDir}/../../sub-${subid}_ses-PNC1_task-idemo_space-MNI152NLin2009cAsym_atlas/sub-${subid}_ses-PNC1_task-idemo_space-MNI152NLin2009cAsym_aal116.nii.gz -expr 'amongst(a,7101,7102)' -prefix ${outDir}/thalamus.nii.gz
# 	3dmaskave -mask $outDir/thalamus.nii.gz -q ${dataFile}> ${outDir}/thalamus_timecourse.1D
	
	3dcalc -a ${dataDir}/../../sub-${subid}_ses-PNC1_task-idemo_space-MNI152NLin2009cAsym_atlas/sub-${subid}_ses-PNC1_task-idemo_space-MNI152NLin2009cAsym_aal116.nii.gz -expr 'amongst(a,7101)' -prefix ${outDir}/L_thalamus.nii.gz 
	3dmaskave -mask $outDir/L_thalamus.nii.gz -q ${dataFile}> ${outDir}/L_thalamus_timecourse.1D
	
	3dcalc -a ${dataDir}/../../sub-${subid}_ses-PNC1_task-idemo_space-MNI152NLin2009cAsym_atlas/sub-${subid}_ses-PNC1_task-idemo_space-MNI152NLin2009cAsym_aal116.nii.gz -expr 'amongst(a,7102)' -prefix ${outDir}/R_thalamus.nii.gz 
	3dmaskave -mask $outDir/R_thalamus.nii.gz -q ${dataFile}> ${outDir}/R_thalamus_timecourse.1D
	
	#Using 3dTcorr1D
	# qsub-run 3dTcorr1D -mask ${dataDir}/../../norm/sub-${subid}_ses-PNC1_task-idemo_space-MNI152NLin2009cAsym_maskStd.nii.gz -prefix ${outDir}/thalamus_conn.nii.gz -overwrite ${dataFile} ${outDir}/thalamus_timecourse.1D
	# 3dcalc -a ${outDir}/thalamus_conn.nii.gz -expr 'atanh(a)' -prefix ${outDir}/thalamus_conn_fisherz.nii.gz -overwrite
	
	qsub-run 3dTcorr1D -mask ${dataDir}/../../norm/sub-${subid}_ses-PNC1_task-idemo_space-MNI152NLin2009cAsym_maskStd.nii.gz -prefix ${outDir}/L_thalamus_conn.nii.gz -overwrite ${dataFile} ${outDir}/L_thalamus_timecourse.1D
	3dcalc -a ${outDir}/L_thalamus_conn.nii.gz -expr 'atanh(a)' -prefix ${outDir}/L_thalamus_conn_fisherz.nii.gz -overwrite
	
	qsub-run 3dTcorr1D -mask ${dataDir}/../../norm/sub-${subid}_ses-PNC1_task-idemo_space-MNI152NLin2009cAsym_maskStd.nii.gz -prefix ${outDir}/R_thalamus_conn.nii.gz -overwrite ${dataFile} ${outDir}/R_thalamus_timecourse.1D
	3dcalc -a ${outDir}/R_thalamus_conn.nii.gz -expr 'atanh(a)' -prefix ${outDir}/R_thalamus_conn_fisherz.nii.gz -overwrite
	

done <$sublist
