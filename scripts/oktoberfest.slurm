#!/bin/bash
#SBATCH -J ok_sb_b1
#SBATCH -o ok_sb_b1_%A_%a.out
#SBATCH -e ok_sb_b1_%A_%a.err
#SBATCH --cpus-per-task=4         		 # should be max(cpus_on_node, number of spectra files)
#SBATCH --partition=shared-cpu
#SBATCH --mem=16G                 		 # ~3GB per spectra file in parallel, i.e. 3*CPUs
#SBATCH --time=4-00:00:00			# should be enough if CPUs = number of spectra files
#SBATCH --get-user-env
#SBATCH --export=NONE
#SBATCH --mail-user=sara.aman@tum.de
#SBATCH --array=0-208%9		 # this email will receive notification about the job
# ================ SETUP ================ #
. /cmnfs/home/s.aman/miniconda3/etc/profile.d/conda.sh
conda activate oktoberfest
module load percolator/3.6.1 thermorawfileparser/1.4.3
# DO NOT CHANGE THIS, only change with --cpus-per-task above!
NUM_THREADS=$SLURM_CPUS_PER_TASK
# =============== CONFIG ================ #
input_dir="/cmnfs/data/proteomics/metaproteomics/PXD011515_sara_analysis/sprotB"
# List of patient IDs
PATIENT_IDS=(11 26 42 49 70 72 77 83 92 93 103 112 120 124 130 135 146 152 153 154 155 161 173 191 194 199 210 222 228 229 231 242 245 251 253 254 256 263 265 266 268 269 276 284 286 295 299 303 310 313 315 318 321 328 334 345 347 351 352 353 360 364 368 373 378 379 380 381 383 400 401 402 416 418 421 422 431 438 442 446 456 459 460 464 465 466 469 470 475 479 480 481 485 491 492 494 497 500 501 503 506 507 508 510 514 515 516 518 520 521 522 524 527 529 530 532 537 538 539 540 542 546 548 550 558 564 565 571 576 578 582 587 588 596 601 604 607 631 637 642 643 648 651 656 662 666 672 680 684 694 696 709 712 718 723 726 730 731 738 741 742 743 754 757 763 764 766 769 770 771 774 780 781 784 785 788 791 802 804 805 810 813 818 824 827 833 839 840 844 848 849 858 859 860 861 866 868 870 874 876 877 887 888 891 895 913 915 923 931)

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
