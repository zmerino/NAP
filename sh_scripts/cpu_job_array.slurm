#!/bin/sh

#SBATCH --partition=Orion
#SBATCH --job-name=cpu_40
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=40        # cpu-cores per task (>1 if multi-threaded tasks)
#SBATCH --mem-per-cpu=4G
#SBATCH --time=48:00:00 

### email notification

#SBATCH --mail-user=zmerino@uncc.edu
#SBATCH --mail-user=zmerino2718@gmail.com
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL

### setup job array
#SBATCH --array=1-10
#SBATCH --output=job_array_40cpu_%j.out

# Specify the path to the config file
config=/users/zmerino/nse_versions/nse_current/sh_scripts/config_archive_40.txt

# Extract the variables for current $SLURM_ARRAY_TASK_ID
distribution_vector=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $2}' $config)
names=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $3}' $config)
trials=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $4}' $config)
min_pow=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $5}' $config)
max_pow=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $6}' $config)
cpu_n=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $7}' $config)
cpp_code=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $8}' $config)

# Print to a file a message that includes the current $SLURM_ARRAY_TASK_ID, the same name, and the sex of the sample
echo "This is array task ${SLURM_ARRAY_TASK_ID}, for variables: ${distribution_vector}, ${names}, ${trials}, ${max_pow}, ${min_pow}, ${cpu_n}, ${cpp_code}" >> job_array_output.txt

# Make temporary directory to store parpool information necessary to run a parallel job.

mkdir -p /users/$USER/.matlab/temp_cluster_jobs/$SLURM_JOB_ID

### select version of Matlab
VER=R2022b

### run job
module load matlab/${VER}
cd /users/zmerino/nse_versions/nse_current

### Change Directory and run MainScript
matlab -nodisplay -nosplash -r "driver_func ${distribution_vector} ${names} ${trials} ${min_pow} ${max_pow} ${cpu_n} ${cpp_code}"

rm -R /users/$USER/.matlab/temp_cluster_jobs/$SLURM_JOB_ID

