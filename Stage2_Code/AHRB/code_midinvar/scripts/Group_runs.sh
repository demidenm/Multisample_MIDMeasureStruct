#!/bin/bash
#
#SBATCH --job-name=AHRB_seclvl_rn
#SBATCH --array=01# run job array
#SBATCH --time=00:20:00
#SBATCH --cpus-per-task=6
#SBATCH --mem-per-cpu=6GB
#SBATCH -p ${PROFILE}
# Outputs ----------------------------------
#SBATCH --output=../logs/Grp_rn.%A_%a.out
#SBATCH --error=../logs/Grp_rn.%A_%a.err
#SBATCH --mail-user=${USER}@stanford.edu
#SBATCH --mail-type=ALL
# ------------------------------------------

source $HOME/miniconda3/etc/profile.d/conda.sh
conda activate fmri_env

# example from job array, sub=("21" "31" "78" "55" "106")
job_id=$( printf %02d ${SLURM_ARRAY_TASK_ID} )
echo "Job ID: " $job_id


# Define directories
dir="/oak/../AHRB/derivatives/analyses/proj_midinvar"
fixed_in="${dir}/firstlvl"
tmp_dir="${dir}/firstlvl_tmp"

# make tmp dir
[ ! -d ${tmp_dir} ] && echo "group out directory exists" | mkdir -p ${tmp_dir}

# filter fixed_in by subj list
sub_list="/oak/../AHRB/code/proj_midinvar/scripts/qc_out/ahrb_final_subjs.tsv"
# cp subjs to tmpdir
cat ${sub_list} | while read line ; do
        sub=$(echo $line | tr -d '"' ) ;
        cp -r ${fixed_in}/${sub} ${tmp_dir}
done

for run in 01 02 ; do

	group_scratch="/scratch/../AHRB/analyses/proj_midinvar/group/run_${run}"
	group_out="${dir}/group/run_${run}"

	[ ! -d ${group_scratch} ] && echo "group scratch directory exists" | mkdir -p ${group_scratch}
	[ ! -d ${group_out} ] && echo "group out directory exists" | mkdir -p ${group_out}

	# run python script
	echo "Subject List: ${sub_list}"
	echo "Run: ${run} "
	echo "Fixed lvl input: ${tmp_dir}"
	echo "Group scratch out: ${group_scratch}"
	echo "Group oak out: ${group_out}"

	python ./python_scripts/mid_group_runs.py ${tmp_dir} ${group_scratch} ${run}

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

done

# remove firstlvl maps in tmp
rm -r ${tmp_dir}