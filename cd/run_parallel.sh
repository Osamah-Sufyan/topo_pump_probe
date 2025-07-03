#!/bin/bash

for param in $(seq 0 -0.025 -0.4); do
    echo "Running with param=$param"
    python pump_probe_cd_t2.py "$param" &
done

wait  # Waits for all background jobs to finish
