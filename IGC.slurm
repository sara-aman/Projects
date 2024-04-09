#!/bin/bash
#SBATCH -J IGC_andromeda_index
#SBATCH -o IGC_andromeda_index_%A_%a.out
#SBATCH -e IGC_andromeda_index_%A_%a.err
#SBATCH --get-user-env
#SBATCH --partition=compms-cpu-big
#SBATCH --nodes=1
#SBATCH --mem=16GB
#SBATCH --cpus-per-task=4
#SBATCH --export=NONE
#SBATCH --mail-user=sara.aman@tum.de
#SBATCH --time=30-00:00:00
#================== PARAMS ===================#
# Change This
MQPAR_DIR="/cmnfs/data/proteomics/metaproteomics/mdpi_proteomes_2019_PXD011515/IGC/4"
#================== SETUP ===================#
#maxquant module load, it will automatically load the required dotnet version
module load maxquant/2.4.4.0
#try this time to load dotnet...
module load dotnet/3.1.426
#================== EXE =====================#
# List of patient IDs
#PATIENT_IDS=(4 9)

#dir="$MQPAR_DIR/${PATIENT_IDS[$SLURM_ARRAY_TASK_ID]}"
#echo "$dir/mqpar.xml"
echo "$MQPAR_DIR/mqpar.xml"

dotnet MaxQuantCmd.exe "$MQPAR_DIR/mqpar.xml"

echo "Done."
