#!/bin/bash

for param in $(seq 0 0.025 0.1); do
    echo "Running with param=$param"
    python current_A0_left_t2_0.py "$param" &
done

wait  # Waits for all background jobs to finish
