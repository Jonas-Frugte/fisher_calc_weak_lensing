#!/bin/bash
#SBATCH --job-name=JFCMBN
#SBATCH --output=/home3/p319950/ResearchProject/outputs/cmb_noise_output.log
#SBATCH --error=/home3/p319950/ResearchProject/outputs/cmb_noise_error.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --time=01:00:00
#SBATCH --mem=10GB

echo "Starting"
source /home3/p319950/venvs/mainenv/bin/activate
cd /home3/p319950/ResearchProject/fisher_calc_weak_lensing/code
python3 -u cmb_noise_fast_runner.py
echo "Done"