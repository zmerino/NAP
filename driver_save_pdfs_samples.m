
% use when driver.m is a script and not a function
clc;clear all; close all;

addpath("functions/")
addpath("compile_nmem_mv/")

% error handling
status = mkdir('log');
diary(fullfile('log','error_log_cpu_failure_distance.txt'))
diary on;

% empty text file used to track progress
filename = ['cpu_failure_distance_script_run-',datestr(datetime(floor(now),'ConvertFrom','datenum')),'.txt'];
full_file = fullfile('log',filename);

fid = fopen(full_file, 'w');
fprintf(fid,['Cpu failure distance script started on: ',datestr(datetime(now,'ConvertFrom','datenum')),'/n']);
fclose(fid);

dir_name = fullfile('data','combined_pdf_estimates');
status = mkdir(dir_name);

% class assignment
actual = distributions;
actual.generate_data = false;

% User Options ============================================================
% script switching board
data_type_flag =            true;   %<- true/false integer powers of 2/real powers of 2
% rndom data generation parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
max_pow =                   20; %<---- maximum exponent to generate samples
min_pow =                   8; %<---- minimum exponent to generate samples
trials =                    30   ;  %<--- trials to run to generate heuristics for programs
step =                      1;  %<---- control synthetic rndom samples to skip being created
temp_min_limit =            0; %<---- set upper limit for both
actual.min_limit =          temp_min_limit;  %<--- lower limit to plot
temp_max_limit =            1000; %<---- set upper limit for both
actual.max_limit =          temp_max_limit; %<--- upper limit to plot
x_resolution =              1000;
% changes how data files are accessed when using data generated from the
% Univariant Random Sample Generator available on zmerino's github
cpu_type =                   '\';%<--- '\' or '/' for windows or linux

distribution_vector = ["Trimodal-Normal","Uniform","Normal","Uniform-Mix","Beta-a0p5-b1p5","Beta-a2-b0p5","Beta-a0p5-b0p5","Generalized-Pareto","Stable"];
distribution = distribution_vector';
names = ["Tri-Modal-Normal","Uniform", "Normal","Uniform-Mix", "Beta(0.5,1.5)", "Beta(2,0.5)", "Beta(0.5,0.5)", "Generalized-Pareto", "Stable"];

% distribution_vector = ["Normal"];
% distribution = distribution_vector';
% names = ["Normal"];


distribution_vector = ["Trimodal-Normal","Uniform","Normal","Uniform-Mix","Beta-a0p5-b1p5","Beta-a2-b0p5","Beta-a0p5-b0p5","Generalized-Pareto"];
distribution = distribution_vector';
names = ["Tri-Modal-Normal","Uniform", "Normal","Uniform-Mix", "Beta(0.5,1.5)", "Beta(2,0.5)", "Beta(0.5,0.5)", "Generalized-Pareto"];


% distribution_vector = ["Beta-a0p5-b1p5","Beta-a2-b0p5","Beta-a0p5-b0p5","Generalized-Pareto"];
% distribution = distribution_vector';
% names = ["Beta(0.5,1.5)", "Beta(2,0.5)", "Beta(0.5,0.5)", "Generalized-Pareto"];
% 
% 
% distribution_vector = ["Beta-a0p5-b1p5"];
% distribution = distribution_vector';
% names = ["Beta(0.5,1.5)"];


% distribution_vector = ["Beta-a0p5-b1p5","Beta-a2-b0p5","Beta-a0p5-b0p5","Generalized-Pareto","Stable"];
% distribution = distribution_vector';
% names = ["Beta(0.5,1.5)", "Beta(2,0.5)", "Beta(0.5,0.5)", "Generalized-Pareto", "Stable"];


% distribution_vector = ["Normal","Uniform"];
% distribution = distribution_vector';
% names = ["Normal","Uniform"];



% find amy of the strings in "str" inside of "distribtuionVector"
str = ["Beta-a0p5-b1p5","Beta-a2-b0p5","Beta-a0p5-b0p5"];
flag = zeros(1,length(distribution_vector));
for a = 1:size(distribution_vector,2)
    for b = 1:size(str,2)
        if strcmp(distribution_vector(a),str(b))
            flag(a) = 1;
            break
        else
            flag(a) = 0;
        end
    end
end

% only works for division with no remainders
fail_nmem = zeros((max_pow-min_pow+1)/step, trials,length(distribution_vector));
fail_nse = zeros((max_pow-min_pow+1)/step, trials,length(distribution_vector));
lagrange_nse = zeros((max_pow-min_pow+1)/step, trials,length(distribution_vector));
lagrange_nmem = zeros((max_pow-min_pow+1)/step, trials,length(distribution_vector));

% table to store meta data
global_table = table();
mse_table = table();

% change code to multithreaded version
for j = 1:length(distribution_vector)
    % Define plot vector for dist_list from 0-1
    if flag(j)
        actual.min_limit = 0;
        actual.max_limit = 1;
        actual.x = linspace(actual.min_limit,actual.max_limit,x_resolution); %<----******** need to update code
    else
        % Define plot vector for distribution from actual.min_limit-actual.max_limit
        actual.min_limit = 0;
        actual.max_limit = 10;
        actual.x = linspace(actual.min_limit,actual.max_limit,x_resolution);
    end
    % Current distribution name
    actual.dist_name = distribution_vector(j);
    % file name for actual distribution. "A_" puts at the top of the folder.
    actual.filename = sprintf(['A_', char(actual.dist_name),'_Act']);
    
    % creat rndom object
    rndom = actual;
    
    T_track = zeros(max_pow-min_pow,1);
    BR_track = zeros(max_pow-min_pow,trials);
    BR0_track = zeros(max_pow-min_pow,trials);
    sample_data = zeros(max_pow-min_pow,1);
    cpu_time_se = [];
    cpu_time_nmem = [];
    sample_track = [];
    
    % Create vector of  samples
    sample_vec = misc_functions.sample_pow(min_pow,max_pow,data_type_flag,step);
    cpu_vec_se = zeros(length(sample_vec),trials);
    cpu_vec_nmem = zeros(length(sample_vec),trials);
    kl_vec_se = zeros(length(sample_vec),trials);
    kl_vec_nmem = zeros(length(sample_vec),trials);
    mse_vec_se = zeros(length(sample_vec),trials);
    mse_vec_nmem = zeros(length(sample_vec),trials);
    block_scale = zeros(4,length(sample_vec),trials);
    block_size = zeros(4,length(sample_vec),trials);

    % store mse per block for all trials per distribution and sample size
    mse_dists_nse = cell(length(sample_vec),trials);
    mse_dists_nmem = cell(length(sample_vec),trials);

    % store block size and scale per distribution and sample size
    block_scale_all = cell(length(sample_vec),trials);
    block_size_all = cell(length(sample_vec),trials);

    actual = actual.dist_list();

    for k = 1:length(sample_vec)

        % generate empty cell arrays to save pdf estiamtes and samples
        
        nse_pdf_data = {};
        nmem_pdf_data = {};

        for i = 1:trials
                    
            % initialize to 0. used to count number of failed NMEM pdfe
            fail_flag = 0;
            
            rndom.Ns = sample_vec(k);
            realx = linspace(actual.min_limit,actual.max_limit,rndom.Ns);
            % p-vector definition for Rtree
            p = [1,0.55,1,0.33,2,ceil(0.0625*rndom.Ns^0.5),40];
            
            % Create rndom.filename for each distribtuion
            rndom.filename = sprintf(['D_', char(actual.dist_name),'_T_','%d', '_S_','%d'],i, rndom.Ns);
            
            %==========================================================
            % % % % % % % start of estimates % % % % % % % % %
            %==========================================================
            tic

            %-- NSE start
            sendFileName1 = ['D_',char(actual.dist_name),cpu_type,char(rndom.filename),'.dat'];
            rndom.randomVSactual = "random";
            rndom = dist_list(rndom);
            sample = rndom.rndData;
            
            tintialSE = cputime;
            % NSE ---------------------------------------------------------
            % nse object instantiation
            nse = NSE;
            nse.max_bs = 1e3;
            nse = nse.stitch(sample);
            
            % extract relevant parameters from object after stich() method
            fail_code = nse.failed;
            SE_x = nse.sx;
            SE_pdf = nse.sPDF;
            SE_cdf = nse.sCDF;
            SE_u = nse.u;
            SE_SQR = nse.sqr;
            nBlocks = nse.nBlocks;
            rndom.Ns = nse.N;
            binrndom.Ns =  nse.binN;
            LG = nse.LG;
            max_LG = nse.LG_max;
            sum_LG = nse.LG_sum;
            T = nse.T;
            BRlevel = nse.BRlevel;
            BR0 = nse.BR0;

            tcpuSE = cputime-tintialSE;

            fail_nse(k,i,j) = fail_code;

            % store information about block
            block_size(1,k,i) = max(nse.block_size);
            block_size(2,k,i) = min(nse.block_size);
            block_size(3,k,i) = mean(nse.block_size);
            block_size(4,k,i) = median(nse.block_size);

            block_scale(1,k,i) = max(nse.block_scale);
            block_scale(2,k,i) = min(nse.block_scale);
            block_scale(3,k,i) = mean(nse.block_scale);
            block_scale(4,k,i) = median(nse.block_scale);

            % NMEM --------------------------------------------------------
            try
                tintialNMEM = cputime;
                [failed, x_NMEM, pdf_NMEM, cdf_NMEM,sqr_NMEM, ~,lagrange_multipler] = EstimatePDF(sample);
                fail_nmem(k,i,j) = 0;
                tcpuNMEM = cputime-tintialNMEM;
            
                n = length(sqr_NMEM);
                dx = 1 / (n + 1);
                u_NMEM = dx:dx:(n * dx);
            catch
                warning('Problem using function.  Assigning a value of 0.');
                lagrange_multipler = 0;
                x_NMEM = linspace(min(sample),max(sample),length(sample))';
                pdf_NMEM = 0*ones(length(sample),1);
                cdf_NMEM = 0*ones(length(sample),1);
                fail_nmem(k,i,j) = 1;
                % don't record cpu time if estimator failed
                tcpuNMEM = NaN;
            end
            
            toc
            %==========================================================
            % % % % % % % % end of estimate % % % % % % % % %
            %==========================================================
            
            %%%%%%%%%%%%%%%% calculate LG multipliers %%%%%%%%%%%%%%%%%%%%%
            lagrange_nse(k,i,j) = sum(max_LG);
            lagrange_nmem(k,i,j) = sum(length(lagrange_multipler));

            %%%%%%%%%%%%%%%%%%%%% store meta data %%%%%%%%%%%%%%%%%%%%%%%%%

            % read in actual distribution
            rndom.dist_list();
            
            % Create final answer file
            rndom.filename = sprintf(['D_', char(rndom.dist_name), ...
                '_T_','%d', '_S_','%d'],i, rndom.Ns);
%             
%             disp([char(actual.dist_name),...
%                 ', Trial: ',num2str(i),'/', num2str(trials), ...
%                 ' sample size: ',num2str(sample_vec(k)), ...
%                 ' failSE: ', num2str(fail_nse(k,i,j)), ' failNMEM: ', ...
%                 num2str(fail_nmem(k,i,j))])

            % store time of computation
            cpu_vec_se(k,i) = tcpuSE;
            cpu_vec_nmem(k,i) = tcpuNMEM;


            
            disp([char(actual.dist_name),...
                ', Trial: ',num2str(i),'/', num2str(trials), ...
                ' sample size: ',num2str(sample_vec(k))])
            

            disp(['estimate size: ',num2str(size(SE_pdf))])

            nse_pdf_data{1,i} = sample;
            nse_pdf_data{2,i} = SE_x;
            nse_pdf_data{3,i} = SE_pdf;

            nmem_pdf_data{1,i} = sample;
            nmem_pdf_data{2,i} = x_NMEM;
            nmem_pdf_data{3,i} = pdf_NMEM;
            
            
        end
        
        filename = [char(actual.dist_name),...
                '_t_', num2str(trials), ...
                '_s_',num2str(sample_vec(k))];

        save(fullfile(dir_name,['nse_', filename, '.mat']), 'nse_pdf_data')
        save(fullfile(dir_name,['nmem_', filename, '.mat']), 'nmem_pdf_data')
    end

    % cpu time
    stdCPUtimeSE = horzcat(sample_vec',std(cpu_vec_se,[],2));
    stdCPUtimeNMEM = horzcat(sample_vec',std(cpu_vec_nmem,[],2));
    avgCPUtimeSE = horzcat(sample_vec',mean(cpu_vec_se,2));
    avgCPUtimeNMEM = horzcat(sample_vec',mean(cpu_vec_nmem,2));
        
    cpu.stdTimeSE{j} = horzcat(sample_vec',std(cpu_vec_se,[],2));
    cpu.stdTimeNMEM{j} = horzcat(sample_vec',std(cpu_vec_nmem,[],2));
    cpu.timeSE{j} = horzcat(sample_vec',mean(cpu_vec_se,2));
    cpu.timeNMEM{j} = horzcat(sample_vec',mean(cpu_vec_nmem,2));

    cpu_nse = cpu_vec_se;
    cpu_nmem = cpu_vec_nmem;

    %%%%%%%%%%%%%%%%%%%%%%% build data table %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % CPU
    temp = vertcat(misc_functions.reshape_groups(sample_vec',cpu_vec_se),...
    misc_functions.reshape_groups(sample_vec',cpu_vec_nmem));

    sample_power = temp(:,1);
    cpu_time = temp(:,2);

    % Distributions 
    distribution = repelem(distribution_vector(j), length(temp(:,2)))';
    name = repelem(names(j), length(temp(:,2)))';

    nse_label = repelem(["NSE"], size(misc_functions.reshape_groups(sample_vec',cpu_vec_se), 1));
    nmem_label = repelem(["NMEM"], size(misc_functions.reshape_groups(sample_vec',cpu_vec_nmem), 1));

    % Failed
    temp = vertcat(misc_functions.reshape_groups(sample_vec',fail_nse(:,:,j)),...
        misc_functions.reshape_groups(sample_vec',fail_nmem(:,:,j)));
 
    fail = temp(:,2);   

    % Lagragian
    temp = vertcat(misc_functions.reshape_groups(sample_vec',lagrange_nse(:,:,j)),...
        misc_functions.reshape_groups(sample_vec',lagrange_nmem(:,:,j)));
    lagrange = temp(:,2);    
    
    estimator = vertcat(nse_label', nmem_label');

    
    % padding NMEM with NaNs because there are no scale/size infromation
    padding = misc_functions.reshape_groups(sample_vec',NaN(length(sample_vec),trials));


    % block size ------------------------

    % max
    bs_max_mat = squeeze(block_size(1,:,:));
    bs_max = misc_functions.reshape_groups(sample_vec',bs_max_mat);
    max_size = vertcat(bs_max, padding);

    sample_power = max_size(:,1);
    max_size = max_size(:,2);


    % min
    bs_min_mat = squeeze(block_size(2,:,:));
    bs_min = misc_functions.reshape_groups(sample_vec',bs_min_mat);
    min_size = vertcat(bs_min, padding);
    min_size = min_size(:,2);


    % mean
    bs_mean_mat = squeeze(block_size(3,:,:));
    bs_mean = misc_functions.reshape_groups(sample_vec',bs_mean_mat);
    mean_size = vertcat(bs_mean, padding);
    mean_size = mean_size(:,2);


    % median
    bs_med_mat = squeeze(block_size(4,:,:));
    bs_med = misc_functions.reshape_groups(sample_vec',bs_med_mat);
    median_size = vertcat(bs_med, padding);
    median_size = median_size(:,2);


    %std dev
    bs_stdev_mat = squeeze(block_size(4,:,:));
    bs_stdev = misc_functions.reshape_groups(sample_vec',bs_stdev_mat);
    std_size = vertcat(bs_stdev, padding);
    std_size = std_size(:,2);

    blocksize = table(max_size, min_size, mean_size, median_size, std_size);


    % block scale ------------------------

    % max
    bs_max_mat = squeeze(block_scale(1,:,:));
    bs_max = misc_functions.reshape_groups(sample_vec',bs_max_mat);
    max_scale = vertcat(bs_max, padding);
    max_scale = max_scale(:,2);

    % min
    bs_min_mat = squeeze(block_scale(2,:,:));
    bs_min = misc_functions.reshape_groups(sample_vec',bs_min_mat);
    min_scale = vertcat(bs_min, padding);
    min_scale = min_scale(:,2);

    % mean
    bs_mean_mat = squeeze(block_scale(3,:,:));
    bs_mean = misc_functions.reshape_groups(sample_vec',bs_mean_mat);
    mean_scale = vertcat(bs_mean, padding);
    mean_scale = mean_scale(:,2);

    % median
    bs_med_mat = squeeze(block_scale(4,:,:));
    bs_med = misc_functions.reshape_groups(sample_vec',bs_med_mat);
    median_scale = vertcat(bs_med, padding);
    median_scale = median_scale(:,2);

    %std dev
    bs_stdev_mat = squeeze(block_size(4,:,:));
    bs_stdev = misc_functions.reshape_groups(sample_vec',bs_stdev_mat);
    std_scale = vertcat(bs_stdev, padding);
    std_scale = std_scale(:,2);

    blockscale = table(max_scale, min_scale, mean_scale, median_scale, std_scale);

    dist_table = table(distribution, name, estimator, sample_power, cpu_time, fail, lagrange, blocksize, blockscale);
    dist_table = splitvars(dist_table);

    % append table per distribution to global table containing data for all
    % distributions
    global_table = [global_table; dist_table];

end

writetable(global_table,fullfile('data','combined_estimator_meta_data.dat'))



