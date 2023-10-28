#!/bin/bash
dir=`pwd`
out_file=${dir}/../../out_midinvar/fixedeff_map_qc.csv

folder=/scratch.global/${USER}/analyses/fixedeff
background_mni=../ROI/MNI152NLin2009cAsym_res-02_desc-brain_T1w.nii.gz
# location nucleus accumbens
vox_x=44
vox_y=72
vox_z=37

min_thresh=0.5
max_thresh=5.0

for sub in `echo ${folder}/sub-* ` ; do
	sub=$(echo $sub | awk -F"/" '{ print $6 }' )
	echo $sub
	
	file=$(echo ${folder}/${sub}/${sub}_ses-2YearFollowUpYArm1_task-MID_effect-fixed_contrast-Lgain-Neut_stat-tstat.nii.gz )
	fsleyes -vl ${vox_x} ${vox_y} ${vox_z} ${background_mni} ${file} -cm red-yellow -dr ${min_thresh} ${max_thresh} 
	echo "Is activation map good (g) or bad (b)?"
	read response_img
	quality=${response_img}
	echo -e "${sub}\t${quality}" >> $out_file
	echo
done
