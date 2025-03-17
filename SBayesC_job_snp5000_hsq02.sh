#!/bin/bash
#SBATCH --job-name=10kh02
#SBATCH --mem=40G
#SBATCH --time=260:30:00
#SBATCH --partition=general
#SBATCH --array=1-100  # Each task corresponds to a file
#SBATCH --output=results_out_snp5k_hsq02/slurm-%j_%a.out  # Save SLURM output in results_out/
#SBATCH --error=results_out_snp5k_hsq02/slurm-%j_%a.err   # Save SLURM errors in results_out/

module load gcc/12.3.0
module load eigen/.3.4.0-gcccore-12.3.0
module load boost/1.82.0-gcc-12.3.0
module load openmpi/4.1.5-gcc-12.3.0
module load r/4.2.1-foss-2022a
export OMP_NUM_THREADS=8

EXECUTABLE="./sbayesc"
EXECUTABLE_NZP="./sbayesc_nzp"
EXECUTABLE_ISS="./sbayesc_ISS"

N_ITER=50000
NN_ITER=50000
PI_INIT=0.01
NUM_SNP=5000
HSQ=0.02
NUM_SEED=1  
R=0.3
NUM_CAUSL=100

TASK_ID=${SLURM_ARRAY_TASK_ID} 
DATA_DIR="data_new/change_hsq/hsq02/"

SIZE_N_SMALL=("500" "5000")

mkdir -p "$DATA_DIR"
for SIZE in "${SIZE_N_SMALL[@]}"; do
    mkdir -p "$DATA_DIR/$SIZE"
done
mkdir -p "$DATA_DIR"
for SIZE in "${SIZE_N_LARGE[@]}"; do
    mkdir -p "$DATA_DIR/$SIZE"
done


for SIZE in "${SIZE_N_SMALL[@]}"; do
    DATA_DIR_READ="data_new/change_hsq/hsq02/$SIZE"
    DATA_DIR_STORE="data_new/change_hsq/hsq02/$SIZE"

    #echo "Generating simulation data..."
    #Rscript generate_simulation.R "$DATA_DIR_STORE" "$HSQ" "$R" "$NUM_SEED" "$NUM_CAUSL" "$SIZE" "$NUM_SNP" "$TASK_ID"
    
    #echo "Running SBayesC_ISS_nu2 on $DATA_DIR (Task ID: $TASK_ID)..."
    #$EXECUTABLE_ISS $DATA_DIR_READ/ldm_data${TASK_ID}.ma \
    #    $DATA_DIR_READ/GWASss_data${TASK_ID}.ma \
    #    $DATA_DIR_STORE/SBayesC_ISS_nu2_result_${TASK_ID}_result.bin \
    #    $DATA_DIR_STORE/SBayesC_ISS_nu2_nnz_ssq_result${TASK_ID}.bin \
    #    $((TASK_ID + 1000)) $SIZE $N_ITER $PI_INIT $HSQ

    #echo "Running SBayesC.R on $DATA_DIR (Task ID: $TASK_ID)..."
    #Rscript sbayesc_r.R $DATA_DIR_READ/ldm_data${TASK_ID}.ma \
    #    $DATA_DIR_READ/GWASss_data${TASK_ID}.ma \
    #    $DATA_DIR_STORE/R_result_${TASK_ID} \
    #    $DATA_DIR_STORE/R_results_nnz_${TASK_ID} \
    #    $((TASK_ID + 1000)) $PI_INIT $HSQ $SIZE $N_ITER

    echo "Running SBayesC on $DATA_DIR (Task ID: $TASK_ID)..."
    $EXECUTABLE $DATA_DIR_READ/ldm_data${TASK_ID}.ma \
        $DATA_DIR_READ/GWASss_data${TASK_ID}.ma \
        $DATA_DIR_STORE/SBayesC_result_${TASK_ID}.bin \
        $DATA_DIR_STORE/SBayesC_nnz_ssq_result${TASK_ID}.bin \
        $((TASK_ID + 1000)) $SIZE $N_ITER $PI_INIT $HSQ

    #echo "Running SBayesC_NZP on $DATA_DIR (Task ID: $TASK_ID)..."
    #$EXECUTABLE_NZP $DATA_DIR_READ/ldm_data${TASK_ID}.ma \
    #    $DATA_DIR_READ/GWASss_data${TASK_ID}.ma \
    #    $DATA_DIR_STORE/SBayesC_ISS_result_${TASK_ID}_result.bin \
    #    $DATA_DIR_STORE/SBayesC_ISS_nnz_ssq_result${TASK_ID}.bin \
    #    $((TASK_ID + 1000)) $SIZE $N_ITER $PI_INIT $HSQ

done


echo "All tasks for changing n completed."


#./sbayesc_ISS data_new/change_hsq/snp10khsq02/500/ldm_data1.ma \
#        data_new/change_hsq/snp10khsq02/500/GWASss_data1.ma \
#        data_new/change_hsq/snp10khsq02/500/SBayesC_ISS_result_1_result.bin \
#        data_new/change_hsq/snp10khsq02/500/SBayesC_ISS_nnz_ssq_result1.bin \
#        $((1 + 1000)) 500 100000 0.02 0.02