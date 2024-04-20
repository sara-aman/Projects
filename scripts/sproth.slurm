#!/bin/bash
#SBATCH -J sproth_b4
#SBATCH -o sproth_b4_%A_%a.out
#SBATCH -e sproth_b4_%A_%a.err
#SBATCH --get-user-env
#SBATCH --partition=shared-cpu
#SBATCH --nodes=1
#SBATCH --mem=16GB
#SBATCH --cpus-per-task=4
#SBATCH --export=NONE
#SBATCH --mail-user=sara.aman@tum.de
#SBATCH --time=4-00:00:00
#SBATCH --array=0-9%5
#================== PARAMS ===================#
# Change This
MQPAR_DIR="/cmnfs/data/proteomics/metaproteomics/PXD011515_sara_analysis/sprotH"
#================== SETUP ===================#
#maxquant module load, it will automatically load the required dotnet version
module load maxquant/2.4.4.0
#try this time to load dotnet...
#module load dotnet/3.1.426
#================== EXE =====================#
# List of patient IDs
PATIENT_IDS=(508 510 514 539 540 542 546 548 870 874) 

dir="$MQPAR_DIR/${PATIENT_IDS[$SLURM_ARRAY_TASK_ID]}"
echo "$dir/mqpar.xml"

dotnet MaxQuantCmd.exe "$dir/mqpar.xml"

echo "Done."
