#!/bin/bash
#SBATCH --job-name=JFFish
#SBATCH --output=/home3/p319950/ResearchProject/outputs/output.log
#SBATCH --error=/home3/p319950/ResearchProject/outputs/error.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --time=24:00:00
#SBATCH --mem=100GB

echo "Starting"
source /home3/p319950/venvs/mainenv/bin/activate
cd /home3/p319950/ResearchProject
python3 -u Fisher_matrix.py # turn cpus per task to 64 and time to 24 hours if this is on
python3 -u Fisher_matrix_powerspectra.py
echo "Done"