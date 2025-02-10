#!/bin/sh
#SBATCH -A b1045
#SBATCH --partition=buyin 
#SBATCH --job-name=cmaq_2016_baseline_validation        # Job name
#SBATCH --output=output_%j.txt       # Output file name (%j will be replaced by the job ID)
#SBATCH --error=error_%j.txt         # Error file name (%j will be replaced by the job ID)
#SBATCH --time=04:00:00              # Time limit (hh:mm:ss)
#SBATCH -n 10                        # Run on a single task (i.e., 1 CPU)
#SBATCH --mail-type=END,FAIL         # Notifications for job done & fail
#SBATCH --mail-user=kyanshlipak2026@u.northwestern.edu  # Your email address for notifications

# Load any necessary modules (adjust to your environment)
#module load python/3.x.x             # Load Python module (version may vary)

# Activate your virtual environment (if needed)
source activate WRF

# Run your Python script
#python cmaq_epa_baseline_kyan_neighbors.py
#python cmaq_epa_2020_RWC.py
#python step1.p
python check_WRF_stats.py
