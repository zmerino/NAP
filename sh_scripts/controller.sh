#!/bin/sh

task1=$(sbatch --parsable --job-name="${JobID1}" cpu_job_array.sh)
sleep 2
echo ${task1}
task1=${task1} 

task2=$(sbatch --parsable --dependency=aftercorr:$task1 --job-name="${JobID2}" post_process.sh)
sleep 2
echo ${task2}