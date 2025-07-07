#!/bin/bash
#SBATCH --job-name=JFFish
#SBATCH --output=/home3/p319950/ResearchProject/fisher_calc_weak_lensing/code/outputs/fish_mat_approx_output.log
#SBATCH --error=/home3/p319950/ResearchProject/fisher_calc_weak_lensing/code/outputs/fish_mat_approx_error.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --time=48:00:00
#SBATCH --mem=70GB

echo "Starting"
source /home3/p319950/venvs/mainenv/bin/activate
cd /home3/p319950/ResearchProject/fisher_calc_weak_lensing/code
python setup.py build_ext --inplace --force
python3 -u Fisher_matrix_approx.py
python3 -u Fisher_matrix_powerspectra.py
echo "Done"