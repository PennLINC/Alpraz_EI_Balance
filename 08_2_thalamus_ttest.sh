#!/bin/bash
nCores=20
sublist="/cbica/projects/alpraz_EI/input/PNC_sublist_FD_age_distance_filter.csv"

scan_file="/cbica/projects/alpraz_EI/input/PNC_scan_list.txt"
echo -n > $scan_file

cov_file="/cbica/projects/alpraz_EI/input/PNC_covariates.csv"
echo "subject FD sex age distance" > $cov_file

mask="/cbica/projects/alpraz_EI/input/PNC_mask_thr995_reslice.nii.gz"

sed 1d $sublist |while IFS=, read subid sesid age sex FD distance other; do
	echo $subid $sesid

	dataDir="/cbica/projects/alpraz_EI/data/PNC/IDEMO_ACOMPCOR_GSR/${subid}/${sesid}"
	
	# Append info to the scan list and covariates files.
	sub_label="sub-${subid}" # This is the label for the scan
	echo "$sub_label ${dataDir}/R_thalamus_conn_fisherz.nii.gz" >> $scan_file
	echo "$sub_label $FD $sex $age $distance" >> $cov_file

done

export OMP_NUM_THREADS=20
# Get the scan list of scans
scans=$(cat $scan_file)

# Call to 3dttest++
# -prefix: output filename
# -overwrite: overwrite the output file if it exists
# -mask: give the mask file: THIS IS IMPORTANT
# -setA label scanlist: we are labeling this "cbf" but can pick anything, also passing our scan list
# -ClustSim #cores: This is a critical step (but takes a while). Can paralellize by setting #cores (20 here)
## Highly recommend reading the help for this fx to understand ClustSim, but in brief:
## ClustSim uses random permutations of your dataset to simulate the null model for your cluster correction.
## It is super convenient because 1) it does everything for you, and 2) it appends the output to your output .nii.gz file which makes it super easy to visualize your results.
echo $nCores
# Change output filename depending on if we are smoothing

3dttest++ -prefix "/cbica/projects/alpraz_EI/output/PNC_ttests/R_thalamus_conn_distance_clustsim.nii.gz" -overwrite -mask $mask -setA thalamusFC $scans -covariates $cov_file -Clustsim 20
