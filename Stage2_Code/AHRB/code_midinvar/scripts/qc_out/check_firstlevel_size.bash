#!/bin/bash

for file_path in $(ls ../../logs/FirstLvl*.err ) ; do 
	file_size=$(stat -c%s "$file_path")
	name=$(basename $file_path) 
	if [ "$file_size" -eq 510 ]; then 
		echo "${name}: exactly 510 bytes." 
	else 
		echo "${name} is NOT 510 bytes. Size: ${file_size}." 
	fi
done
