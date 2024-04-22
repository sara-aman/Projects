#!/bin/bash
#SBATCH -J ok_sb_766
#SBATCH -o ok_sb_766_%A_%a.out
#SBATCH -e ok_sb_766_%A_%a.err
#SBATCH --cpus-per-task=4         		 # should be max(cpus_on_node, number of spectra files)
#SBATCH --partition=shared-cpu
#SBATCH --mem=16G                 		 # ~3GB per spectra file in parallel, i.e. 3*CPUs
#SBATCH --time=4-00:00:00			# should be enough if CPUs = number of spectra files
#SBATCH --get-user-env
#SBATCH --export=NONE
#SBATCH --mail-user=sara.aman@tum.de
#SBATCH --array=0		 # this email will receive notification about the job
# ================ SETUP ================ #
. /cmnfs/home/s.aman/miniconda3/etc/profile.d/conda.sh
conda activate oktoberfest
module load percolator/3.6.1 thermorawfileparser/1.4.3
# DO NOT CHANGE THIS, only change with --cpus-per-task above!
NUM_THREADS=$SLURM_CPUS_PER_TASK
# =============== CONFIG ================ #
input_dir="/cmnfs/data/proteomics/metaproteomics/PXD011515_sara_analysis/sprotB"
# List of patient IDs
PATIENT_IDS=(766)

cat > rescoring_config.json << EOL
{
    "type": "Rescoring",
    "tag": "",
    "output": "$input_dir/${PATIENT_IDS[$SLURM_ARRAY_TASK_ID]}/oktoberfest_output",
    "inputs": {
        "search_results": "$input_dir/${PATIENT_IDS[$SLURM_ARRAY_TASK_ID]}/combined/txt/msms.txt",
        "search_results_type": "Maxquant",
        "spectra":"$input_dir/${PATIENT_IDS[$SLURM_ARRAY_TASK_ID]}",
        "spectra_type": "raw"
    },
    "models": {
        "intensity": "Prosit_2020_intensity_HCD",
        "irt": "Prosit_2019_irt"
    },
    "prediction_server": "koina.wilhelmlab.org:443",
    "numThreads": $NUM_THREADS,
    "fdr_estimation_method": "percolator",
    "allFeatures": false,
    "regressionMethod": "spline",
    "ssl": true,
    "thermoExe": "/cmnfs/home/s.aman/ThermoRawFileParser1.4.3/ThermoRawFileParser.exe",
    "massTolerance": 20,
    "unitMassTolerance": "ppm",
    "ce_alignment_options": {
        "ce_range": [19,50],
        "use_ransac_model": false
    }
}
EOL

# =========== Run oktoberfest =========== #

/usr/bin/time -v -o stats.txt /cmnfs/home/s.aman/miniconda3/envs/oktoberfest/bin/python -m oktoberfest --config_path=rescoring_config.json

exit 0
