#!/bin/bash
#
#SBATCH --job-name=MLS_firstlvl
#SBATCH --array=1,3,8,12-18,22,26,28,30-34,37,39-41,46,48,50,51,53,56,58,60-63,67,72,75-79,81-83,85,87,89-91,93-95,98-99,101-102,104-105,109-110,112-114,116,119,121-122,127-132,136,139,144,149,151,155-156,158-160,162,164-165,167-168,175-176,183,185,195,197,199-200,204,207,213,220,224,226-227,230,238%20
#SBATCH --time=00:12:00
#SBATCH --cpus-per-task=6
#SBATCH --mem-per-cpu=4GB
#SBATCH -p ${PROFILE}
# Outputs ----------------------------------
#SBATCH --output=../logs/FirstLvl.%A_%a.out
#SBATCH --error=../logs/FirstLvl.%A_%a.err
#SBATCH --mail-user=${USER}@stanford.edu
#SBATCH --mail-type=ALL
# ------------------------------------------

source $HOME/miniconda3/etc/profile.d/conda.sh
conda activate fmri_env

# Define directories
data_in="/oak/../data/MLS"
data_preproc="${data_in}/derivatives/fmriprep_23.1.0"

# example from job array, sub=("21" "31" "78" "55" "106")
subj=$( printf %02d ${SLURM_ARRAY_TASK_ID} )
echo "SUBJECT_ID: " $subj
sub="sub-${subj}"

scratch_out="/scratch/../MLS/analyses/proj_midinvar/firstlvl/${sub}"
analysis_out="${data_in}/derivatives/analyses/proj_midinvar/firstlvl/${sub}"

[ ! -d ${scratch_out} ] && echo "scratch directory exists" | mkdir -p ${scratch_out}
[ ! -d ${analysis_out} ] && echo "scratch directory exists" | mkdir -p ${analysis_out} 

echo "#### Starting First level run GLMs ####"
# run python script
echo "sub: ${sub}"
echo "data_in: ${data_in}"
echo "data_preproc: ${data_preproc}"
echo "scratch_out: ${scratch_out}"
echo "analysis_out: ${analysis_out}"


python ./mid_firstlevel.py ${sub} ${data_in} ${data_preproc} ${scratch_out} ${analysis_out}

firstlvl_error=$?

if [ ${firstlvl_error} -eq 0 ]; then
        echo "Python first level completed successfully!"
else
    	echo "Python first level failed."
        exit 1
fi

echo "1. Syncing files from scratch to analysis path. Deleted from scratch once sync'd"
 
rsync -av --remove-source-files ${scratch_out}/ ${analysis_out}/

echo
echo
echo "#### Starting Precision Weighted Fixed Effect ####" 

# Define directories
firstlvl_inp="${data_in}/derivatives/analyses/proj_midinvar/firstlvl/${sub}"
fixed_scratch="/scratch/../MLS/analyses/proj_midinvar/fixedeff/${sub}"
fixed_out="${data_in}/derivatives/analyses/proj_midinvar/fixedeff/${sub}"

[ ! -d ${fixed_scratch} ] && echo "scratch directory exists" | mkdir -p ${fixed_scratch}
[ ! -d ${fixed_out} ] && echo "fixed out directory exists" | mkdir -p ${fixed_out} 

# run python script
echo "sub: ${sub}"
echo "First lvl input: ${firstlvl_inp}"
echo "Fixed scratch out: ${fixed_scratch}"
echo "Fixed oak out: ${fixed_out}"

python ./python_scripts/mid_fixedeff.py ${sub} ${firstlvl_inp} ${fixed_scratch} 

fixeff_error=$?
if [ ${fixeff_error} -eq 0 ]; then
        echo "Python fixed effect completed successfully!"
else
    	echo "Python first effect failed."
        exit 1
fi

echo "Syncing files from scratch to analysis path. Deleted from scratch once sync'd"
 
rsync -av --remove-source-files ${fixed_scratch}/ ${fixed_out}/

echo "2. Syncing files from scratch to analysis path. Deleted from scratch once sync'd"
