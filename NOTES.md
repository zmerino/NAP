

Raw data location on the cluster ---------------------------

/users/zmerino/nse_versions/nse_current/


Raw data directories ---------------------------------------

data_cpu_20/

coontains an old run comparing NAPS and NMEM.
100 trials.
12 cpus used.

data_cpu_20_wall/

coontains an old run comparing NAPS and NMEM. 
100 trials.
12 cpus used.
Wall time data included.

data_cpu_40_wall/

coontains an old run comparing NAPS and NMEM. 
100 trials.
40 cpus used.
Wall time data included.

data_large_n/
1 trials.
40 cpus used.
Wall time data included.

Things to keep in mind -----------------------------------

1) each subdirectory in the parent data directory contains
    raw estimates per sample size for all trials.
    
2) each .dat/.mat file in the parent directory contains
    meta data for each distribution i.e. cpu times, 
    lagrange multipliers, sample size, etc.

3) to process the meta data run the process_meta_data.ma
    script and point to the desired directory.

4) to process the estimate data for KL/MSE/SQR/PDFS/etc
    run the process_kl_mse_quantiles.m script and point
    to a directory that contains all of the estiamte data
    files. NOTE: you must copy the data in each 
    subdirectory it a single directory.

5) NOTE: the process_kl_mse_quantiles.m script is quite
    slow for large sample sizes and even slower if
    computing the MSE/KL per quantile.
