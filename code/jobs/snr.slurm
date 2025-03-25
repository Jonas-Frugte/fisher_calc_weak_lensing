#!/bin/bash
#SBATCH --job-name=JFSNR
#SBATCH --output=/home3/p319950/ResearchProject/outputs/snr_output.log
#SBATCH --error=/home3/p319950/ResearchProject/outputs/snr_error.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --time=6:00:00
#SBATCH --mem=200GB
# used to be 64 cpus per task and 6 hours
echo "Starting"
source /home3/p319950/venvs/mainenv/bin/activate
cd /home3/p319950/ResearchProject
python3 -u SNR.py
echo "Done"