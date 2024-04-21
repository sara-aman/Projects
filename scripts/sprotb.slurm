#!/bin/bash
#SBATCH -J sprotb_b5
#SBATCH -o sprotb_b5_%A_%a.out
#SBATCH -e sprotb_b5_%A_%a.err
#SBATCH --get-user-env
#SBATCH --partition=shared-cpu
#SBATCH --nodes=1
#SBATCH --mem=16GB
#SBATCH --cpus-per-task=4
#SBATCH --export=NONE
#SBATCH --mail-user=sara.aman@tum.de
#SBATCH --time=4-00:00:00
#SBATCH --array=0-77%10
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
PATIENT_IDS=(4 9 42 49 72 93 103 242 263 276 286 295 299 303 310 313 315 318 321 328 334 345 347 351 352 353 368 373 401 416 431 465 466 470 475 480 485 491 492 494 497 501 516 521 522 527 530 558 565 566 601 637 642 643 712 731 764 766 769 770 774 780 784 785 802 824 833 839 840 844 849 859 860 866 895 913 923 931)


dir="$MQPAR_DIR/${PATIENT_IDS[$SLURM_ARRAY_TASK_ID]}"
echo "$dir/mqpar.xml"

dotnet MaxQuantCmd.exe "$dir/mqpar.xml"

echo "Done."
