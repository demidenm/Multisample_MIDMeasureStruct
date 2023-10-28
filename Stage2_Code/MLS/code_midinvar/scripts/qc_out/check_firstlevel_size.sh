#!/bin/bash

for file_path in $(ls ../../logs/FirstLvl*.err ) ; do 
	file_size=$(stat -c%s "$file_path")
	name=$(basename $file_path) 
	if [ "$file_size" -ne 510 ] && [ "$file_size" -ne 1124 ] ; then 
		echo "${name} is NOT 510 or 1124 bytes. Size: ${file_size}." 
	fi
done
