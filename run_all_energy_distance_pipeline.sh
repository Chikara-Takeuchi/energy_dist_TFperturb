#!/usr/bin/bash
#SBATCH -J edist_pipeline_step1_2      # Job name
#SBATCH -N 1                          # Total number of nodes requested (16 cores/node)
#SBATCH -t 24:00:00                   # Run time (hh:mm:ss) - 20 hrs limit
#SBATCH -p GPUv100s
#SBATCH -o run_output.out
#SBATCH -e run_output.err
#SBATCH --gres=gpu:1

TARGET_FILE_ID="syn70753570"
SYNAPSE_TOKEN=${1}

CONTAINER_PATH="./edist_pipeline.sif"
CONFIG_PATH="./data/config1_2.json"
BIN_PATH="./energy_dist_pipeline/bin"

#If data folder doesn't exist, prepare data folder.
mkdir -p "./data"

apptainer pull ${CONTAINER_PATH} docker://docker.io/takechikara/energy_distance_env:latest
git clone https://github.com/Chikara-Takeuchi/energy_dist_pipeline.git

echo "Preprocess mudata file"
apptainer exec --nv ${CONTAINER_PATH} python preprocess_mudata.py --synapse_id ${TARGET_FILE_ID} --auth_token ${SYNAPSE_TOKEN}

echo "[Step1] Filtering outlier gRNAs"
apptainer exec --nv ${CONTAINER_PATH} python ${BIN_PATH}/1_filtering_gRNA.py ${CONFIG_PATH}

echo "[Step2] calculate energy distance between targets and non-targeting"
apptainer exec --nv ${CONTAINER_PATH} python ${BIN_PATH}/2_e_distance_nontargeting.py ${CONFIG_PATH}

echo "[Step2_1] visualize results of energy distance analysis"
apptainer exec --nv ${CONTAINER_PATH} python ${BIN_PATH}/2_1_Plot_figure.py ${CONFIG_PATH}

echo "All steps completed successfully."