#!/bin/bash
#SBATCH --job-name=JFDatExp
#SBATCH --output=/home3/p319950/ResearchProject/fisher_calc_weak_lensing/code/outputs/data_export_output.log
#SBATCH --error=/home3/p319950/ResearchProject/fisher_calc_weak_lensing/code/outputs/data_export_error.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --time=6:00:00
#SBATCH --mem=60GB

echo "Starting"
source /home3/p319950/venvs/mainenv/bin/activate
cd /home3/p319950/ResearchProject/fisher_calc_weak_lensing/code
export OMP_NUM_THREADS=1
export OMP_DYNAMIC=FALSE
python3 -u data_exporter_new.py