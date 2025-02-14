#!/bin/bash

# Load required modules
module load eigen/.3.4.0-gcccore-12.3.0
module load r/4.2.1-foss-2022a

# Define paths
DATA_DIR="data/large_samples/1e7"
mkdir -p "$DATA_DIR"  # Create the directory if it doesn't exist

# Step 1: Generate 100 datasets using R
echo "Generating datasets using R..."
Rscript generate_simulation.R "$DATA_DIR"

# Step 2: Define executable and parameters
EXECUTABLE="./sbayesc"
N_ITER=10000
PI_INIT=0.1
HSQ_INIT=0.5

# Step 3: Submit jobs to run SBayesC for each dataset
qsubshcom "$EXECUTABLE $DATA_DIR/ldm_data{TASK_ID}.ma $DATA_DIR/GWASss_data{TASK_ID}.ma \
    $DATA_DIR/ldm_data{TASK_ID}_result.bin $DATA_DIR/nnz_ssq_result{TASK_ID}.bin $((1000 + {TASK_ID})) \
    $N_ITER $PI_INIT $HSQ_INIT" \
    8 8G sbayesc_task 0:30:00 "-array=1-100"

echo "All jobs submitted!"
