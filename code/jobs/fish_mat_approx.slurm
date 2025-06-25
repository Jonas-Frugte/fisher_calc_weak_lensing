#!/bin/bash
#SBATCH --job-name=JFFish
#SBATCH --output=/home3/p319950/ResearchProject/fisher_calc_weak_lensing/code/outputs/fish_mat_approx_output_exp_fix.log
#SBATCH --error=/home3/p319950/ResearchProject/fisher_calc_weak_lensing/code/outputs/fish_mat_approx_error_exp_fix.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --time=5-00:00:00
#SBATCH --mem=70GB

echo "Starting"
source /home3/p319950/venvs/mainenv/bin/activate
cd /home3/p319950/ResearchProject/fisher_calc_weak_lensing/code
python3 -u Fisher_matrix_approx.py # turn cpus per task to 64 and time to 24 hours if this is on
#python3 -u Fisher_matrix_powerspectra.py
echo "Done"