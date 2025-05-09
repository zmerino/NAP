

% this script is used for locally testing scripts used on the cluster.
% note that lines 

% define variables that the slurm script calls from the config.txt file
% in the folder sh_scripts.

distribution_vector = ["Beta-a0p5-b0p5"];
names = ["Beta(0.5,0.5)"];
distribution_vector = ["Normal"];
names = ["Normal"];
trials = '3';
min_pow = '11';
max_pow = '11';
cpu_n = '20';
% cpp_code = "cpp_code/";
cpp_code = "cpp_code_smooth/";


% CALL FUCNTIONS THAT COMPARES NAPS (serial/parallel) AND NMEM

% driver_func(distribution_vector,names,trials, min_pow,max_pow,cpu_n,cpp_code)

% CALL FUNCTION THAT RUNS NAPS FOR LARGE SAMPLE SIZE

% driver_func(distribution_vector,names,trials, min_pow,max_pow,cpu_n,cpp_code)
driver_func_n(distribution_vector,names,trials, min_pow,max_pow,cpu_n,cpp_code)
% driver_func_kl(distribution_vector,names,trials, min_pow,max_pow,cpu_n,cpp_code)

% read_path = "data_cpu_20_wall_2";
% % process_func(distribution_vector,names,trials, min_pow,max_pow,cpu_n,read_path)
% process_func_nse(distribution_vector,names,trials, min_pow,max_pow,cpu_n,read_path)
% process_func_nmem(distribution_vector,names,trials, min_pow,max_pow,cpu_n,read_path)

% CALLS FUNCTION FOR MULIPLE DISTRIBUTIONS
% mimics how slurm submits a job array
% NOT REALLY NEEDED, ABOVE FUNCTIONS ARE FINE FOR TESTING ON LOCAL PC

% distribution_vector = ["Normal","Beta-a0p5-b0p5"];
% names = ["Normal", "Beta(0.5,0.5)"];
% 
% trials = '2';
% min_pow = '13';
% max_pow = '14';
% cpu_n = '20';
% cpp_code = "cpp_code/";
% 
% for i =1:length(distribution_vector)
%     driver_func(distribution_vector(i),names(i),trials, min_pow,max_pow,cpu_n,cpp_code)
% end