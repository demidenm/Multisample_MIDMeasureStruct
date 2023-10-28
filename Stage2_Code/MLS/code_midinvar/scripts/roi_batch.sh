#!/bin/bash
#
#SBATCH --job-name=MLS_ROIest
#SBATCH --array=1 #Only needs 1 job, calculates ROI for all maps within input dir 
#SBATCH --time=00:25:00
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=8GB
#SBATCH -p russpold,normal,owners
# Outputs ----------------------------------
#SBATCH --output=../logs/ROI.%A_%a.out
#SBATCH --error=../logs/ROI.%A_%a.err
#SBATCH --mail-user=demidenm@stanford.edu
#SBATCH --mail-type=ALL
# ------------------------------------------

source $HOME/miniconda3/etc/profile.d/conda.sh
conda activate fmri_env

# Define directories
data_in="/oak/stanford/groups/russpold/data/MLS"
roi_dir="${data_in}/code/proj_midinvar/ROI"
first_in="${data_in}/derivatives/analyses/proj_midinvar/firstlvl"
roi_scratch="/scratch/groups/russpold/MLS/analyses/proj_midinvar/ROI_est"
roi_out="${data_in}/derivatives/analyses/proj_midinvar/ROI_est"
data_preproc="${data_in}/derivatives/fmriprep_23.1.0"

# example from job array, sub=("21" "31" "78" "55" "106")
#subj=$( printf %02d ${SLURM_ARRAY_TASK_ID} )
echo "Starting job."


[ ! -d ${roi_scratch} ] && echo "scratch ROI dir exists" | mkdir -p ${roi_scratch}
[ ! -d ${roi_out} ] && echo "ROI directory exists" | mkdir -p ${roi_out} 

echo "#### Starting Script to Extract ROIs ####"
# run python script
echo "data_in: ${data_in}"
echo "roi_dir: ${roi_dir}" 
echo "first_in: ${first_in}"
echo "roi_scratch: ${roi_scratch}"  
echo "roi_out: ${roi_out}"
echo "data_preproc: ${data_preproc}"


python ./betamap_roi.py ${data_preproc} ${roi_dir} ${first_in} ${roi_scratch}

roi_error=$?
if [ ${roi_error} -eq 0 ]; then
        echo "Python ROI script completed successfully!"
else
    	echo "Python ROI script failed."
        exit 1
fi


echo "Syncing files from scratch to analysis path. Deleted from scratch once sync'd"
 
rsync -av --remove-source-files ${roi_scratch}/ ${roi_out}/

echo "DONE"
