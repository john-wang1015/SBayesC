#!/bin/bash
#SBATCH --job-name=sbayesc_array
#SBATCH --mem=20G
#SBATCH --time=30:30:00
#SBATCH --partition=general
#SBATCH --array=1-100  # Each task corresponds to a file
#SBATCH --output=results_out/slurm-%j_%a.out  # Save SLURM output in results_out/
#SBATCH --error=results_out/slurm-%j_%a.err   # Save SLURM errors in results_out/


module load gcc/12.3.0
module load eigen/.3.4.0-gcccore-12.3.0
module load r/4.2.1-foss-2022a

# Executables
EXECUTABLE="./sbayesc"
EXECUTABLE_NZP="./sbayesc_nzp"
N_ITER=10000
PI_INIT=0.1
HSQ_INIT=0.5

TASK_ID=${SLURM_ARRAY_TASK_ID}  # Slurm assigns this automatically

# List of dataset sizes
DATA_SIZES=("1e4" "1e5" "1e6" "1e7" "1e8" "1e9")

# Ensure required directories exist before running the R script
for SIZE in "${DATA_SIZES[@]}"; do
    if [[ "$SIZE" =~ ^1e[789]$ ]]; then
        mkdir -p "data/large_samples/$SIZE"
    else
        mkdir -p "data/samll_samples/$SIZE"
    fi
done

# Generate simulation data using R script
echo "Generating simulation data..."
Rscript generate_simulations.R "${DATA_SIZES[@]}"

# Process each dataset
for SIZE in "${DATA_SIZES[@]}"; do
    if [[ "$SIZE" =~ ^1e[789]$ ]]; then
        DATA_DIR="data/large_samples/$SIZE"
    else
        DATA_DIR="data/samll_samples/$SIZE"
    fi

    # Run SBayesC.R
    echo "Running SBayesC.R on $DATA_DIR (Task ID: $TASK_ID)..."
    Rscript SBayesC.R $DATA_DIR/ldm_data${TASK_ID}.ma \
        $DATA_DIR/GWASss_data${TASK_ID}.ma \
        $DATA_DIR/beta_r_results_${TASK_ID}.ma \
        $DATA_DIR/nnz_r_results_${TASK_ID}.ma \
        $((TASK_ID + 1000)) $SIZE $N_ITER

    echo "Running SBayesC on $DATA_DIR (Task ID: $TASK_ID)..."
    $EXECUTABLE $DATA_DIR/ldm_data${TASK_ID}.ma \
        $DATA_DIR/GWASss_data${TASK_ID}.ma \
        $DATA_DIR/ldm_data${TASK_ID}_result.bin \
        $DATA_DIR/nnz_ssq_result${TASK_ID}.bin \
        $((TASK_ID + 1000)) $N_ITER $PI_INIT $HSQ_INIT

    echo "Running SBayesC_NZP on $DATA_DIR (Task ID: $TASK_ID)..."
    $EXECUTABLE_NZP $DATA_DIR/ldm_data${TASK_ID}.ma \
        $DATA_DIR/GWASss_data${TASK_ID}.ma \
        $DATA_DIR/ldm_data${TASK_ID}_nzp_result.bin \
        $DATA_DIR/nnz_ssq_result${TASK_ID}_nzp.bin \
        $((TASK_ID + 1000)) $SIZE $N_ITER $PI_INIT $HSQ_INIT
done

echo "All tasks completed."
