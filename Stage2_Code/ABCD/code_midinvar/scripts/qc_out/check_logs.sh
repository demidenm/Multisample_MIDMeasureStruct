#!/bin/bash 
for logs in $(echo ../batch_logs/*.out ) ; do
	failed_result=$(grep -r -1 "failed" $logs )
	if [[ -n $failed_result ]]; then
		base_log=$(basename $logs )
    		echo "failed: $base_log "
	fi
done
