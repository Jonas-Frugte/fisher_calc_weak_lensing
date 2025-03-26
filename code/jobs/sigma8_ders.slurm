#!/bin/bash
#SBATCH --job-name=JFS8DER
#SBATCH --output=/home3/p319950/ResearchProject/outputs/sigma8_ders_output.log
#SBATCH --error=/home3/p319950/ResearchProject/outputs/sigma8_ders_error.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=20:00:00
#SBATCH --mem=50GB

echo "Starting"
source /home3/p319950/venvs/mainenv/bin/activate
cd /home3/p319950/ResearchProject/fisher_calc_weak_lensing/code
python3 -u sigma8_ders.py
echo "Done"