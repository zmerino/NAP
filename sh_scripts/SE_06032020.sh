#! /bin/sh

#PBS -N SE_04032020

#PBS -l nodes=1:ppn=20
#PBS -l walltime=300:00:00
#PBS -q copperhead

#PBS -m bea

#PBS -M zmerino@uncc.edu

cd /users/zmerino/SE_06032020

### load Matlab 2018Rb
module load matlab/R2018b

### Change Directory and run MainScript
matlab -nodisplay -nosplash < driver.m > run.log

### exit matlab
#exit