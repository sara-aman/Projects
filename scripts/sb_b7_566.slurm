#!/bin/bash
#SBATCH -J sb_b7_566
#SBATCH -o sb_b7_566_%A.out
#SBATCH -e sb_b7_566_%A.err
#SBATCH --get-user-env
#SBATCH --partition=shared-cpu
#SBATCH --nodes=1
#SBATCH --mem=16GB
#SBATCH --cpus-per-task=4
#SBATCH --export=NONE
#SBATCH --mail-user=sara.aman@tum.de
#SBATCH --time=4-00:00:00
#================== PARAMS ===================#
# Change This
MQPAR_DIR="/cmnfs/data/proteomics/metaproteomics/PXD011515_sara_analysis/sprotB"
#================== SETUP ===================#
#maxquant module load, it will automatically load the required dotnet version
module load maxquant/2.4.4.0
#try this time to load dotnet...
#module load dotnet/3.1.426
#================== EXE =====================#
# List of patient IDs
PATIENT_IDS=(566)

dir="$MQPAR_DIR/$PATIENT_IDS"
echo "$dir/mqpar.xml"

dotnet MaxQuantCmd.exe "$dir/mqpar.xml"

echo "Done."
