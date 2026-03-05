#!/bin/bash
#SBATCH --job-name=JFFish
#SBATCH --output=/home3/p319950/ResearchProject/fisher_calc_weak_lensing/code/outputs/fish_mat_output.log
#SBATCH --error=/home3/p319950/ResearchProject/fisher_calc_weak_lensing/code/outputs/fish_mat_error.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --time=4-00:00:00
#SBATCH --mem=100GB

echo "Starting"
source /home3/p319950/venvs/mainenv/bin/activate
cd /home3/p319950/ResearchProject/fisher_calc_weak_lensing/code
python setup.py build_ext --inplace --force
python3 -u Fisher_matrix.py
echo "Done"