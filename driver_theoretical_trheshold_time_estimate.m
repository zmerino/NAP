% use when driver.m is a script and not a function
clc;clear;close all;

% User Options ============================================================
data_type_flag =            true;   %<- true/false integer powers of 2/real powers of 2

% rndom data generation parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
max_pow =                   26; %<---- maximum exponent to generate samples
min_pow =                   8; %<---- minimum exponent to generate samples
trials =                    1000;  %<--- trials to run to generate heuristics for programs
max_trials = 100;
max_sample = 2^(22);
step =                      1;  %<---- control synthetic rndom samples to skip being created
temp_min_limit =            0; %<---- set upper limit for both
actual.min_limit =          temp_min_limit;  %<--- lower limit to plot
temp_max_limit =            1000; %<---- set upper limit for both
actual.max_limit =          temp_max_limit; %<--- upper limit to plot
x_resolution =              1000;
% changes how data files are accessed when using data generated from the
% Univariant Random Sample Generator available on zmerino's github
cpu_type =                   '\';%<--- '\' or '/' for windows or linux

distribution_vector = ["Uniform-Mix","Generalized-Pareto","Stable","Trimodal-Normal","Normal", "Uniform","Beta-a0p5-b1p5","Beta-a2-b0p5","Beta-a0p5-b0p5"];

% empty text file used to track progress
status = mkdir('log');
filename = ['Threshold_Total_Time-', datestr(datetime(floor(now), 'ConvertFrom', 'datenum')), '.txt'];
full_file = fullfile('log', filename);

fid = fopen(full_file, 'w');
fprintf(fid,['Time estimate for threshold run started on: ', datestr(datetime(now, 'ConvertFrom','datenum')), '\n']);
fprintf(fid,['Sample Range: ', num2str(2^(min_pow)),' - ', num2str(2^(max_pow)), '\n']);
fprintf(fid,['Number of Trials: ', num2str(trials), '\n']);
fclose(fid);

% class assignment
actual = distributions;
actual.generate_data = false;
% NOTE: true_quantile is commented out because it is slow with Stable dist

%% Main Function Call Loop used to lable plot figures

% find any of the strings in "str" inside of "distribtuionVector"
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

% function estimate_time(distribution_vector,trials,min_pow,max_pow)

% Initialize 3D matirx with zeros that is trials by # of samples sizes by
% # number of distributions.
dist_t_track = zeros(trials,(max_pow-min_pow+1),length(distribution_vector));
dist_br0_track = zeros(trials,(max_pow-min_pow+1),length(distribution_vector));

TotalIterations = length(distribution_vector) * trials * (max_pow - min_pow + 1);

total_time_per_dist_track = zeros(1,max_pow - min_pow + 1);
total_time = zeros(1,length(distribution_vector));

tTotalStart = tic;
currentIteration = 1;
% Loop over every distribution type
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
    
    % Initialize empty or zero filled matrices/arrays ---------------------
    
    T_track = NaN(1,max_pow-min_pow+1);
    BR0_track = NaN(1,max_pow-min_pow+1);
    dist_BR0_track = cell(length(distribution_vector));
    dist_BR_track = cell(length(distribution_vector));
    
    % Create vector of  samples
    sample_vec = misc_functions.sample_pow(min_pow,max_pow,data_type_flag,step);
    BoxCPUtimeSE = zeros((max_pow-min_pow)*trials,2);
    BoxCPUtimeNMEM = zeros((max_pow-min_pow)*trials,2);
    
    actual = actual.dist_list();
    
    for k = 1:length(sample_vec)
        
        % initialize to 0. used to count number of failed NMEM pdfe
        fail_flag = 0;
        
        rndom.Ns = sample_vec(k);
        realx = linspace(actual.min_limit,actual.max_limit,rndom.Ns);
        % p-vector definition for Rtree
        p = [1,0.55,1,0.33,2,ceil(0.0625*rndom.Ns^0.5),40];
        
        % Create rndom.filename for each distribtuion
        rndom.filename = sprintf(['D_', char(actual.dist_name),'_T_','%d', '_S_','%d'],i, rndom.Ns);
        
        rndom.randomVSactual = "random";
        rndom = dist_list(rndom);
        sample = rndom.rndData;
        
        tstart = tic;
        
        x = sort(sample);
        
        % Track T,BR per trial
        [pList,T,BRlevel,BR0] = block_definition.r_tree(x,actual.min_limit,actual.max_limit,rndom.filename,p);
        
        tend = toc(tstart);

        PercentComplete = 100 * (currentIteration / TotalIterations);
        currentIteration = currentIteration + 1;
        
        % Estimate total time per distribution and sample for all
        % trials
        h_conversion = 3600;

        % Account for only 100 trials above a certainty sample size,
        % becuase the BR0 converges to zero
        if max_sample <= sample_vec(k)
            total_time_per_sample = tend*max_trials/h_conversion;

            total_time_per_dist_track(1,k) = total_time_per_sample;
    
            meta_string = ['\n', char(actual.dist_name),...
                ' Sample size: ', num2str(sample_vec(k)),...
                ' Total Time per Sample (100 trials): ', num2str(total_time_per_sample), ' hrs'];
        else
            total_time_per_sample = tend*trials/h_conversion;

            total_time_per_dist_track(1,k) = total_time_per_sample;
    
            meta_string = ['\n', char(actual.dist_name),...
                ' Sample size: ', num2str(sample_vec(k)),...
                ' Total Time per Sample: ', num2str(total_time_per_sample), ' hrs'];
        end
        
    
        disp(meta_string)
    
        fid = fopen(full_file, 'at');
        fprintf(fid, meta_string);
        fclose(fid);
            
    end

    total_time_per_dist = sum(total_time_per_dist_track,2);
    total_time(j) = total_time_per_dist;

    meta_string = ['\n', char(actual.dist_name),...
        ' Total Time per Distribution: ', num2str(total_time_per_dist), ' hrs'];
    
    disp(meta_string)

    fid = fopen(full_file, 'at');
    fprintf(fid, meta_string);
    fprintf(fid,'\n_____________________________________________________');
    fclose(fid);
        
end

    meta_string = ['\nTotal Time of Computation: ', num2str(sum(total_time,2)), ' hrs'];
    
    disp(meta_string)

    fid = fopen(full_file, 'at');
    fprintf(fid,'\n');
    fprintf(fid,'\n=====================================================');
    fprintf(fid, meta_string);
    fclose(fid);
