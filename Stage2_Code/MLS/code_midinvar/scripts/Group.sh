#!/bin/bash
#
#SBATCH --job-name=MLS_seclvl
#SBATCH --array=01# run job array
#SBATCH --time=00:20:00
#SBATCH --cpus-per-task=6
#SBATCH --mem-per-cpu=4GB
#SBATCH -p russpold,normal,owners
# Outputs ----------------------------------
#SBATCH --output=../logs/Grp.%A_%a.out
#SBATCH --error=../logs/Grp.%A_%a.err
#SBATCH --mail-user=demidenm@stanford.edu
#SBATCH --mail-type=ALL
# ------------------------------------------

source $HOME/miniconda3/etc/profile.d/conda.sh
conda activate fmri_env

# example from job array, sub=("21" "31" "78" "55" "106")
job_id=$( printf %02d ${SLURM_ARRAY_TASK_ID} )
echo "Job ID: " $job_id


# Define directories
dir="/oak/stanford/groups/russpold/data/MLS/derivatives/analyses/proj_midinvar"
fixed_in="${dir}/fixedeff"
tmp_dir="${dir}/fixedeff_tmp"
group_scratch="/scratch/groups/russpold/MLS/analyses/proj_midinvar/group"
group_out="${dir}/group"

# filter fixed_in by subj list
sub_list="/oak/stanford/groups/russpold/data/MLS/code/proj_midinvar/scripts/qc_out/mls_final_subjs.tsv"

[ ! -d ${group_scratch} ] && echo "group scratch directory exists" | mkdir -p ${group_scratch}
[ ! -d ${group_out} ] && echo "group out directory exists" | mkdir -p ${group_out}
[ ! -d ${tmp_dir} ] && echo "group out directory exists" | mkdir -p ${tmp_dir}

# cp subjs to tmpdir
cat ${sub_list} | while read line ; do 
	sub=$(echo $line | tr -d '"' ) ; 
	cp -r ${fixed_in}/${sub} ${tmp_dir}
done 

# run python script
echo "Subject List: ${sub_list}"
echo "Fixed lvl input: ${tmp_dir}"
echo "Group scratch out: ${group_scratch}"
echo "Group oak out: ${group_out}"

python ./mid_group.py ${tmp_dir} ${group_out}

grp_error=$?
if [ ${grp_error} -eq 0 ]; then
	echo "Python group level completed successfully!"
else
	echo "Python group level failed."
	exit 1
fi

# t to cohens d conversion of maps
python ./python_scripts/tstat_to_cohensd.py ${group_scratch}

echo "Syncing files from scratch to analysis path. Deleted from scratch once sync'd" 
rsync -av --remove-source-files ${group_scratch}/ ${group_out}/

rm -r ${tmp_dir}
