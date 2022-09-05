#!/bin/sh

#SBATCH --partition=Orion
#SBATCH --job-name=nse_threshold
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=48:30:00 

### email notification

#SBATCH --mail-user=zmerino@uncc.edu
#SBATCH --mail-user=zmerino2718@gmail.com
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL

### select version of Matlab

#VER=R2018b
VER=R2018b
#VER=R2020b

### run job

module load matlab/${VER}
cd /users/zmerino/nse_versions/nse_current

### Change Directory and run MainScript

matlab -nodisplay -nosplash < driver_extrema.m > run.log

