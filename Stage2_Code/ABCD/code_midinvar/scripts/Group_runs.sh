#!/bin/bash -l
#SBATCH -J grp_runs
#SBATCH --array=1 # jobs
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=10G
#SBATCH -t 02:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=${USER}.edu
#SBATCH -p msismall,amdsmall
#SBATCH -o batch_logs/%x_%A_%a.out
#SBATCH -e batch_logs/%x_%A_%a.err
#SBATCH -A ${PROFILE}

source $HOME/miniconda3/etc/profile.d/conda.sh
conda activate fmri_env
module load fsl

# Define directories
fixed_in=/scratch.global/{USER}/analyses/firstlvl
group_dir=/scratch.global/{USER}/analyses/group
script_dir="/home/../${USER}/analyses/code_midinvar/scripts"
sample=first1k

for run in 01 02 ; do
	tmp_dir=/scratch.global/{USER}/tmp/fixedeff_${run}
	group_scratch=/tmp/group/run_${run}
	group_out=${group_dir}/run_${run}
	sub_list="${script_dir}/qc_out/subsamp_ids/subs-${sample}.tsv"

	[ ! -d ${group_scratch} ] && echo "group scratch directory exists" | mkdir -p ${group_scratch}
	[ ! -d ${group_out} ] && echo "group out directory exists" | mkdir -p ${group_out}
	[ ! -d ${tmp_dir} ] && echo "group out directory exists" | mkdir -p ${tmp_dir}

	# cp subjs to tmpdir
	cat ${sub_list} | while read line ; do 
		cp -r ${fixed_in}/sub-${line} ${tmp_dir}
	done 

	# run python script
	echo "Sample Type: ${sample}"
	echo "Run: ${run}"
	echo "Subject List: ${sub_list}"
	echo "Fixed lvl input: ${tmp_dir}"
	echo "Group scratch out: ${group_scratch}"
	echo "Group oak out: ${group_out}"

	python ${script_dir}/python_scripts/mid_group_runs.py ${tmp_dir} ${group_out} ${run}
	grp_error=$?

        if [ ${grp_error} -eq 0 ]; then
                echo "Python group level completed successfully!"
        else
    	echo "Python group level failed."
                exit 1
        fi

	echo "Syncing files from scratch to analysis path. Deleted from scratch once sync'd" 
	rsync -av --remove-source-files ${group_scratch}/ ${group_out}/

	rm -r ${tmp_dir}

	# convert run t-maps to cohens d
	echo
	echo
	echo "Converting t-stats to cohen's D effect size maps"
	python ${script_dir}/python_scripts/tstat_to_cohensd.py ${group_out}
done
