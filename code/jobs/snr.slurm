#!/bin/bash
#SBATCH --job-name=JFSNR
#SBATCH --output=/home3/p319950/ResearchProject/fisher_calc_weak_lensing/code/outputs/snr_output.log
#SBATCH --error=/home3/p319950/ResearchProject/fisher_calc_weak_lensing/code/outputs/snr_error.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --time=3:00:00
#SBATCH --mem=200GB
# used to be 64 cpus per task and 6 hours
echo "Starting"
source /home3/p319950/venvs/mainenv/bin/activate
cd /home3/p319950/ResearchProject/fisher_calc_weak_lensing/code
python3 -u SNR.py
echo "Done"