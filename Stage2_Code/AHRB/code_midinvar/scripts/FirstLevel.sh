#!/bin/bash
#
#SBATCH --job-name=AHRB_firstlvl
#SBATCH --array=1-108%20 # sub-26 doesnt have run2
#SBATCH --time=00:20:00
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
data_in="/../AHRB"
data_preproc="${data_in}/derivatives/fmriprep_23.1.0"

# example from job array, sub=("21" "31" "78" "55" "106")
subj=$( printf %02d ${SLURM_ARRAY_TASK_ID} )
echo "SUBJECT_ID: " $subj
sub="sub-${subj}"

scratch_out="/scratch/../AHRB/analyses/proj_midinvar/firstlvl/${sub}"
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


python ./python_scripts/mid_firstlevel.py ${sub} ${data_in} ${data_preproc} ${scratch_out} ${analysis_out}

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
fixed_scratch="/scratch/../AHRB/analyses/proj_midinvar/fixedeff/${sub}"
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
