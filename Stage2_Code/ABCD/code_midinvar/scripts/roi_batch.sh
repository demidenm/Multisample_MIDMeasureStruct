#!/bin/bash

#SBATCH -J roi_est
#SBATCH --array=1 # jobs
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=6
#SBATCH --mem-per-cpu=8G
#SBATCH -t 03:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=${USER}@umn.edu
#SBATCH -p msismall,amdsmall
#SBATCH -o batch_logs/%x_%A_%a.out
#SBATCH -e batch_logs/%x_%A_%a.err
#SBATCH -A ${PROFILE}

source $HOME/miniconda3/etc/profile.d/conda.sh
conda activate fmri_env
module load fsl

# define variables
home_out=/home/faird/${USER}/analyses/out_midinvar
ses=ses-2YearFollowUpYArm1
s3_bucket=s3://ABCD_BIDS/derivatives/fmriprep_v23_1_0
roi_dir=/home/faird/${USER}/analyses/code_midinvar/ROI
tmp_in=/tmp/processed
roi_scratch=/tmp/roi
roi_out=${home_out}/roi
firstlvl_inp=/scratch.global/${USER}/analyses/firstlvl
script_dir=/home/faird/${USER}/analyses/code_midinvar/scripts

# create dirs 
[ ! -d ${roi_scratch} ] && echo "scratch ROI dir exists" | mkdir -p ${roi_scratch}
[ ! -d ${roi_out} ] && echo "ROI directory exists" | mkdir -p ${roi_out}

# copy data from beh ouput & s3, ignore unnecessary files

echo "Copying MID files from s3 bucket"
for subs in $(echo ${firstlvl_inp}/sub-* ) ; do
	sub=$(basename $subs ) 
	mkdir -p ${tmp_in}/${sub}/${ses}/func/

	# get brain mask for calculating overlay %	

	for filename in $(s3cmd ls ${s3_bucket}/${ses}/${sub}/${ses}/func/ \
        	        | grep "MID" \
			| grep "brain_mask.nii.gz" \
			| grep -vE "fsLR|goodvoxels|fsnative" \
			| awk '{ print $NF }') ; do
		s3cmd get ${filename} ${tmp_in}/${sub}/${ses}/func --skip-existing
	done
done

echo "#### Starting Script to Extract ROIs ####"
# run python script
echo "preproc: tmp_in: ${tmp_in}"
echo "roi_dir: ${roi_dir}"
echo "firstlvl_inp: ${firstlvl_inp}"
echo "roi_scratch: ${roi_scratch}"
echo "roi_out: ${roi_out}"


python ${script_dir}/python_scripts/betamap_roi.py ${tmp_in} ${roi_dir} ${firstlvl_inp} ${roi_scratch}
roi_error=$?

if [ ${roi_error} -eq 0 ]; then
        echo "Python ROI script completed successfully!"
else
        echo "Python ROI script failed."
        exit 1
fi
echo
echo "Syncing files from scratch to analysis path. Deleted from scratch once sync'd"
rsync -av --remove-source-files ${roi_scratch}/ ${roi_out}/

