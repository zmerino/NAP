
function process_func_nmem(distribution_vector,names,trials, min_pow,max_pow,cpu_n,read_path)

addpath("functions/")
% addpath("functions_plotting/")

% JOB ARRAY VARIABLES -----------------------------------------------------

% cast characters to strings or to ints for relevant variables
distribution_vector = string(distribution_vector)
names = string(names)
trials = str2num(trials)
min_pow = str2num(min_pow)
max_pow = str2num(max_pow)
cpu_n = str2num(cpu_n)
read_path = string(read_path)

%%%%%%%%%%%%%%%%% COMMENT OUT IF RUNNING LOCALLY %%%%%%%%%%%%%%%%%%%%%%%%%%

% create parallel process on cluster

% create a local cluster object
% pc = parcluster('local')

% explicitly set the JobStorageLocation to the temp directory that was
% created in your sbatch script
% pc.JobStorageLocation = strcat('/home/zmerino/.matlab/temp_cluster_jobs/', getenv('SLURM_JOB_ID'))

% start the matlabpool with maximum available workers
% control how many workers by setting ntasks in your sbatch script
poolobj = gcp('nocreate'); % get the pool object, and do it avoiding creating a new one.
if isempty(poolobj) % check if there is not a pool.
    % create parallel process with given parameters
%     parpool(pc, str2num(getenv('SLURM_CPUS_ON_NODE')), 'IdleTimeout', Inf)
    parpool('Processes', str2num(getenv('SLURM_CPUS_ON_NODE')), 'IdleTimeout', Inf)
else
    delete( gcp('nocreate')); % delete the current pool object.
%     % create parallel process with given parameters
    parpool(pc, str2num(getenv('SLURM_CPUS_ON_NODE')), 'IdleTimeout', Inf)
    parpool('Processes', str2num(getenv('SLURM_CPUS_ON_NODE')), 'IdleTimeout', Inf)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%% To Plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

plot_mse_dist = false;

plot_pdf = true;
plot_cdf = false;
plot_sqr = true;
plot_heavy = false;
save_figs = false;

calc_q = true;

% choose to visualise figures or not
fig_plot = 'on';

%%%%%%%%%%%%%%%%%%%%%%%%%%% Import data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Trials per sample
trials = 100;
plt_trials = 10;
% number of quantile segments
qnum = 50;
% number of cpus used to compute raw data
% cpu_n = 20;

% Create actual distrobution onjects
actual = distributions;

% Path to directory
% dir_name = fullfile('data','pdf_estimates');
table_name = 'mse_kl_100.dat';
% dir_name = fullfile(strcat('data_cpu_',num2str(cpu_n),'_wall_2'));
dir_name = fullfile(read_path);
write_dir = fullfile(dir_name,strcat('kl_mse_cpu_',num2str(cpu_n),'_t_',num2str(trials)));
status = mkdir(write_dir);

% Define the etimates to plot

% Sample range
n_vec = 2.^[min_pow:max_pow];

% d_vec = ["Trimodal-Normal","Uniform","Normal","Beta-a0p5-b1p5","Beta-a2-b0p5","Beta-a0p5-b0p5","Generalized-Pareto","Stable"];
% names = ["Trimodal-Normal","Uniform","Normal","Beta(0.5,1.5)","Beta(2,0.5)", "Beta(0.5,0.5)", "Generalized-Pareto","Stable"];

d_vec = distribution_vector;
% names = names;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% quantiles %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

actual = distributions;
actual.generate_data = false;

% Initialize empty structure to store data in
n = length(n_vec);
d = length(d_vec);
t = length(d_vec);
nap = cell(d,n,t,7);
nmem = cell(d,n,t,7);

mse_dist_nse_per_q = zeros(qnum,plt_trials);
mse_dist_nmem_per_q = zeros(qnum,plt_trials);
kl_dist_nmem_per_q = zeros(qnum,plt_trials);
kl_dist_nse_per_q = zeros(qnum,plt_trials);

global_table = table();

for i = 1:length(d_vec)


    % save distribution name to object
    dist_name = d_vec(i);

    for j = 1:length(n_vec)
        parfor k = 1:trials

            actual = distributions;
            actual.generate_data = false;

            % save distribution name to object
            actual.dist_name = dist_name;

            % Make sure that distribution names are saved as character vectors

            % track calculations
            disp(['dist: ',num2str(i),'/',num2str(length(d_vec)),...
                ' s: ',num2str(j),'/',num2str(length(n_vec)),...
                ' t: ',num2str(k),'/',num2str(trials)])

            % Load data from estimate directory
            read_dir = fullfile(dir_name, strcat(d_vec(i),'_cpu_',num2str(cpu_n),'_t_',num2str(trials)));
            nmem = load(fullfile(read_dir,['nmem_',num2str(d_vec(i)),'_t_',num2str(trials),'_s_', num2str(n_vec(j)),'.mat']));

            % x
            nmem_data = nmem.nmem_pdf_data;

            nmem_x = nmem_data{2,k};

            % pdf
            nmem_pdf = nmem_data{3,k};

            % cdf
            nmem_cdf = utils_analysis.trapz(nmem_data{2,k}, nmem_data{3,k});

            % calcualte kl,mse,quantile -----------------------------------
            % NMEM --------------------------------------------------------
            xs = nmem_x;
            fs = nmem_pdf;
            Fs = nmem_cdf;
            % get actual distribution
            actual.min_limit = min(xs);
            actual.max_limit = max(xs);
            actual.x = xs;
            actual = actual.dist_list();
            %             disp(actual)

            if sum(~isfinite(xs)) || sum(~isfinite(fs))|| sum(~isfinite(Fs))
                warning('non-finite values')
                test = 'test';
            end

            % MSE per quantile
            if calc_q
                [mse_dist_nmem_per_q(:,k), quantiles,xq,fq] = utils_analysis.quant_metric(xs, fs, Fs, qnum, d_vec(i),'MSE');
            end
            % MSE per trial
            mse_dist_nmem(j, k) = utils_analysis.mse(actual.pdf_y, fs);
            % KL per quantile
            if calc_q
                [kl_dist_nmem_per_q(:,k), ~,~,~] = utils_analysis.quant_metric(xs, fs, Fs, qnum, d_vec(i),'KL');
            end
            % KL per trial
            kl_dist_nmem(j, k) = utils_analysis.kl(actual.pdf_y, fs');

        end

        if calc_q
            filename = sprintf('kl_per_q_%s_n_%s_t_%s.mat',d_vec(i), num2str(n_vec(j)), num2str(trials));
            save(fullfile(write_dir,['nmem_', filename]), 'kl_dist_nmem_per_q','-v7.3')

            filename = sprintf('mse_per_q_%s_n_%s_t_%s.mat',d_vec(i), num2str(n_vec(j)), num2str(trials));
            save(fullfile(write_dir,['nmem_', filename]), 'mse_dist_nmem_per_q','-v7.3')
        end


    end

    % build table


    % MSE
    temp = utils.reshape_groups(n_vec',mse_dist_nmem);

    sample_power = temp(:,1);
    mse_total = temp(:,2);

    % KL
    temp = utils.reshape_groups(n_vec',kl_dist_nmem);

    kl_total = temp(:,2);

    % Distributions
    distribution = repelem(d_vec(i), length(temp(:,2)))';
    name = repelem(names(i), length(temp(:,2)))';

    nmem_label = repelem(["NMEM"], size(utils.reshape_groups(n_vec',mse_dist_nmem), 1));


    estimator = nmem_label';


    dist_table = table(distribution, name, estimator, sample_power, mse_total,kl_total);
    dist_table = splitvars(dist_table);

    % append table per distribution to global table containing data
    % for all distributions
    global_table = [global_table; dist_table];
end


writetable(global_table,fullfile(write_dir,table_name),"WriteMode","append")
end