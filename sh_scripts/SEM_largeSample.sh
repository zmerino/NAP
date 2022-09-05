#! /bin/sh

#PBS -N SEM_largeSample

#PBS -l nodes=1:ppn=20
#PBS -l walltime=400:00:00
#PBS -q copperhead


#PBS -m bea

#PBS -M zmerino@uncc.edu

cd /users/zmerino/SEM_largeSample

### load Matlab 2018Rb
module load matlab/R2018b

### Change Directory and run MainScript
matlab -nodisplay -nosplash < stitchPDF_Driver.m > run.log

### exit matlab
#exit