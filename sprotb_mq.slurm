#!/bin/bash
#SBATCH -J sprotb_mq
#SBATCH -o sprotb_mq_%A_%a.out
#SBATCH -e sprotb_mq_%A_%a.err
#SBATCH --get-user-env
#SBATCH --partition=shared-cpu
#SBATCH --nodelist=ernie.exbio.wzw.tum.de
#SBATCH --nodes=1
#SBATCH --mem=16GB
#SBATCH --cpus-per-task=4
#SBATCH --export=NONE
#SBATCH --mail-user=sara.aman@tum.de
#SBATCH --time=4-00:00:00
#SBATCH --array=0-211%5
#================== PARAMS ===================#
# Change This
MQPAR_DIR="/cmnfs/data/proteomics/metaproteomics/mdpi_proteomes_2019_PXD011515/sprotB"
#================== SETUP ===================#
#maxquant module load, it will automatically load the required dotnet version
module load maxquant/2.4.4.0
#try this time to load dotnet...
module load dotnet/3.1.426
#================== EXE =====================#
# List of patient IDs
PATIENT_IDS=(4 9 11 26 42 49 70 72 77 83 92 93 103 112 120 124 130 135 146 152 153 154 155 161 173 191 194 199 210 222 228 229 231 242 245 251 253 254 256 263 265 266 268 269 276 284 286 295 299 303 310 313 315 318 321 328 334 345 347 351 352 353 360 364 368 373 378 379 380 381 383 400 401 402 416 418 421 422 431 438 442 446 456 459 460 464 465 466 469 470 475 479 480 481 485 491 492 494 497 500 501 503 506 507 508 510 514 515 516 518 520 521 522 524 527 529 530 532 537 538 539 540 542 546 548 550 558 564 565 566 571 576 578 582 587 588 596 601 604 607 631 637 642 643 648 651 656 662 666 672 680 684 694 696 709 712 718 723 726 730 731 738 741 742 743 754 757 763 764 766 769 770 771 774 780 781 784 785 788 791 802 804 805 810 813 818 824 827 833 839 840 844 848 849 858 859 860 861 866 868 870 874 876 877 887 888 891 895 913 915 923 931)

dir="$MQPAR_DIR/${PATIENT_IDS[$SLURM_ARRAY_TASK_ID]}"
echo "$dir/mqpar.xml"

dotnet MaxQuantCmd.exe "$dir/mqpar.xml"

echo "Done."
