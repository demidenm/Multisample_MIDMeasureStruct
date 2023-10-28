#!/bin/bash -l

cd $HOME

if [ -d ${HOME}/miniconda3 ] ; then

	echo "miniconda installed in $HOME "
	echo "	... continuing "
else

	echo "miniconda not found in $HOME "
	echo "	... installing "
	
	# download latest miniconda3 env to home directory
	wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh

	# install miniconda
	bash Miniconda3-latest-Linux-x86_64.sh
fi


#if grep -q 'export PATH="$HOME/miniconda3/bin:$PATH"' ${HOME}/.bashrc; then
#
#	echo "miniconda path line already exists in .bashrc."
#else
#
#	echo 'export PATH="$HOME/miniconda3/bin:$PATH"' >> ${HOME}/.bashrc
#	echo "Miniconda kine has been added to .bashrc."
#	rm ${HOME}/Miniconda3-latest-Linux-x86_64.sh
#fi

if [ -d ${HOME}/miniconda3/envs/fmri_env ] ; then

	echo "fmri_env exists within miniconda in home. "

else 
	echo "fmri_env doesnt exist within home, creating and isntall packages."
	# create environment for analysis
	conda create --name fmri_env python=3.9.7
	#conda activate fmri_env

	conda install --name fmri_env --file ${HOME}/analyses/config/analysis_packages.txt

fi
