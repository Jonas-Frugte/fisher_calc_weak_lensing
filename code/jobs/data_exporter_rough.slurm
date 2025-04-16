#!/bin/bash
#SBATCH --job-name=JFDatExp
#SBATCH --output=/home3/p319950/ResearchProject/fisher_calc_weak_lensing/code/outputs/data_export_rough_output.log
#SBATCH --error=/home3/p319950/ResearchProject/fisher_calc_weak_lensing/code/outputs/data_export_rough_error.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=40
#SBATCH --time=5:00:00
#SBATCH --mem=100GB

echo "Starting"

source /home3/p319950/venvs/mainenv/bin/activate

cd /home3/p319950/ResearchProject/fisher_calc_weak_lensing/code
python3 -u data_exporter.py

echo "Done"