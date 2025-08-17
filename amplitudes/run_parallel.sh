#!/bin/bash

for param in $(seq 0 0.025 0.1); do
    echo "Running with param=$param"
    python pump_probe_cd_ampli_topo_left.py "$param" &
done

wait  # Waits for all background jobs to finish
