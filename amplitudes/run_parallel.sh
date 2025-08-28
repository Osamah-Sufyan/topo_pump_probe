#!/bin/bash

LOGFILE="all_output.log"
: > "$LOGFILE"

max_jobs=10
count=0

for param in $(seq 0 0.05 0.4); do
    echo "Running with param=$param" | tee -a "$LOGFILE"
    python "pump_probe_ampli_topo_left.py" "$param" >> "$LOGFILE" 2>&1 &
    ((count++))

    # If we've launched max_jobs, wait for them to finish
    if (( count % max_jobs == 0 )); then
        wait
    fi
done

wait  # wait for remaining jobs
