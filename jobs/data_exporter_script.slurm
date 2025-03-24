#!/bin/bash
#SBATCH --job-name=JFDatExp
#SBATCH --output=/home3/p319950/ResearchProject/outputs/data_export_output.log
#SBATCH --error=/home3/p319950/ResearchProject/outputs/data_export_error.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=27
#SBATCH --time=20:00:00
#SBATCH --mem=160GB

echo "Starting"

source /home3/p319950/venvs/mainenv/bin/activate

python3 -u /home3/p319950/ResearchProject/data_exporter.py

echo "Done"