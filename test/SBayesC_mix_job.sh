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
EXECUTABLE_I="./sbayesc_I"
N_ITER=10000
PI_INIT=0.1
HSQ_INIT=0.5
numSNP = 1000

TASK_ID=${SLURM_ARRAY_TASK_ID}  # Slurm assigns this automatically

# List of dataset sizes
DATA_SIZES=("1e5" "1e6")

# Ensure required directories exist before running the R script
for SIZE in "${DATA_SIZES[@]}"; do
    mkdir -p "data_mix/$SIZE"
done

# Generate simulation data using R script
echo "Generating simulation data..."
Rscript generate_mix_simulation.R "${DATA_SIZES[@]}"

# Process each dataset
for SIZE in "${DATA_SIZES[@]}"; do
    DATA_DIR="data_mix/$SIZE"

    # Run SBayesC.R
    echo "Running SBayesC.R on $DATA_DIR (Task ID: $TASK_ID)..."
    Rscript sbayesc_r.R $DATA_DIR/ldm1_data_mix${TASK_ID}.ma \
        $DATA_DIR/GWASss_data3_mix${TASK_ID}.ma \
        $DATA_DIR/beta_r_results_${TASK_ID} \
        $DATA_DIR/nnz_r_results_${TASK_ID} \
        $((TASK_ID + 1000)) $SIZE $N_ITER

    echo "Running SBayesC on $DATA_DIR (Task ID: $TASK_ID)..."
    $EXECUTABLE $DATA_DIR/ldm1_data_mix${TASK_ID}.ma \
        $DATA_DIR/GWASss_data3_mix${TASK_ID}.ma \
        $DATA_DIR/ldm_data${TASK_ID}_result.bin \
        $DATA_DIR/nnz_ssq_result${TASK_ID}.bin \
        $((TASK_ID + 1000)) $N_ITER $PI_INIT $HSQ_INIT

    echo "Running SBayesC_NZP on $DATA_DIR (Task ID: $TASK_ID)..."
    $EXECUTABLE_NZP $DATA_DIR/ldm1_data_mix${TASK_ID}.ma \
        $DATA_DIR/GWASss_data3_mix${TASK_ID}.ma \
        $DATA_DIR/ldm_data${TASK_ID}_nzp_result.bin \
        $DATA_DIR/nnz_ssq_result${TASK_ID}_nzp.bin \
        $((TASK_ID + 1000)) $SIZE $N_ITER $PI_INIT $HSQ_INIT

    echo "Running SBayesC_I on $DATA_DIR (Task ID: $TASK_ID)..."
    $EXECUTABLE $DATA_DIR/GWASss_data3_mix${TASK_ID}.ma \
        $DATA_DIR/ldm_data${TASK_ID}_result.bin \
        $DATA_DIR/nnz_ssq_result${TASK_ID}.bin \
        $((TASK_ID + 1000)) $N_ITER $PI_INIT $HSQ_INIT $numSNP

done

echo "starting plot....."
Rscript plot_box.R

echo "All tasks completed."

