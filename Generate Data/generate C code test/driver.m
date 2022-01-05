% function driver(percSample,numSubs)
% ^^ for use with Master_driver.m to different subsampling parameters
profile on
% track command window
diary commandWindowOUT.txt
% use when driver.m is a script and not a function
clc;clear all; close all;
% class assignment
actual = distributions;
actual.generate_data = false;

% SUBSAMPLING PARAMETERS---------------------------------------------------

BootSample = "off";
if BootSample == "off"
    % percentage of sample used to create subsample
    percSample = 1;
    % number of subsamples to generate
    numSubs = 1;
else
    % percentage of sample used to create subsample
    percSample = 0.4;
    % number of subsamples to generate
    numSubs = 30;
end

%--------------------------------------------------------------------------

tic
% User Options ============================================================
% script switching board
estimator_call_flag =       true;   %<- true/false call SE on/off
estimator_plot_flag =       false;   %<- true/false plot SE results on/off
data_type_flag =            true;   %<- true/false integer powers of 2/real powers of 2
save_graphics =             false;   %<- true/false save .png of plots on/off
% rndom data generation parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
max_pow =                   11; %<---- maximum exponent to generate samples
min_pow =                   9; %<---- minimum exponent to generate samples
trials =                    3;  %<--- trials to run to generate heuristics for programs
step =                      1;  %<---- control synthetic rndom samples to skip being created
temp_min_limit =            0; %<---- set upper limit for both
actual.min_limit =          temp_min_limit;  %<--- lower limit to plot
temp_max_limit =            1000; %<---- set upper limit for both
actual.max_limit =          temp_max_limit; %<--- upper limit to plot
x_resolution =              1000;
% changes how data files are accessed when using data generated from the
% Univariant Random Sample Generator available on zmerino's github
cpu_type =                   '\';%<--- '\' or '/' for windows or linux

% Example distribution to test %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
distribution_vector = ["Beta-a0p5-b1p5","Beta-a2-b0p5","Beta-a0p5-b0p5",...
    "Bimodal-Normal","BirnbaumSaunders","Burr",...
    "Exponential","Extreme-Value","Gamma","Generalized-Extreme-Value",...
    "Generalized-Pareto","HalfNormal","Normal","Square-periodic",...
    "tLocationScale","Uniform","Uniform-Mix","Weibull","Chisquare",...
    "InverseGaussian","Trimodal-Normal","Stable",...
    "Stable2","Stable3","Stable1","BirnbaumSaunders-Stable"];
% distribution_vector = ["Normal","Uniform","Beta-a0p5-b0p5","Generalized-Pareto","Stable"];
distribution_vector = ["Beta-a0p5-b0p5","Generalized-Pareto"];
% distribution_vector = ["Beta-a0p5-b0p5","Generalized-Pareto";
distribution_vector = ["Trimodal-Normal","Normal","Uniform","Beta-a0p5-b1p5","Beta-a2-b0p5","Beta-a0p5-b0p5","Generalized-Pareto"];
% distribution_vector = ["Beta-a0p5-b0p5"];

distribution_vector = ["Uniform-Mix","Generalized-Pareto","Stable","Trimodal-Normal","Normal", "Uniform","Beta-a0p5-b1p5","Beta-a2-b0p5","Beta-a0p5-b0p5"];
distribution_vector = ["Generalized-Pareto","Stable","Normal", "Uniform","Beta-a0p5-b1p5","Beta-a2-b0p5","Beta-a0p5-b0p5"];
distribution_vector = ["Uniform-Mix","Generalized-Pareto","Stable","Trimodal-Normal","Normal", "Uniform","Beta-a0p5-b1p5","Beta-a2-b0p5","Beta-a0p5-b0p5"];
distribution_vector = ["Trimodal-Normal","Normal", "Uniform"];
% distribution_vector = ["Beta-a0p5-b1p5"];
% distribution_vector = ["Stable"];
% distribution_vector = ["Normal"];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Main Function Call Loop used to lable plot figures

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
fail_nmem = zeros(trials,(max_pow-min_pow)/step,length(distribution_vector));
fail_se = zeros(trials,(max_pow-min_pow)/step,length(distribution_vector));

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
    sample_data_box = [];
    kl_se = [];
    kl_nmem = [];
    Hellinger_se = [];
    Hellinger_nmem = [];
    max_LG_se = [];
    max_LG_nmem = [];
    MSE_se = [];
    MSE_nmem = [];
    all_se = {};
    all_nmem = {};
    cpu_time_se = [];
    cpu_time_nmem = [];
    sample_track = [];
    
    % Create vector of  samples
    sample_vec = misc_functions.sample_pow(min_pow,max_pow,data_type_flag,step);
    cpu_vec_se = zeros(length(sample_vec),trials);
    cpu_vec_nmem = zeros(length(sample_vec),trials);
    BoxCPUtimeSE = [];
    BoxCPUtimeNMEM = [];
    
    for i = 1:trials
        
        actual = actual.dist_list();
        
        for k = 1:length(sample_vec)
            
            interp_etimate = [];
            estimate_data = {};
            
            % initialize to 0. used to count number of failed NMEM pdfe
            fail_flag = 0;
            
            rndom.Ns = sample_vec(k);
            realx = linspace(actual.min_limit,actual.max_limit,rndom.Ns);
            % p-vector definition for Rtree
            p = [1,0.55,1,0.33,2,ceil(0.0625*rndom.Ns^0.5),40];
            % only effects threshold
            p = [1,1,0,0.5,0.5,ceil(0.0625*rndom.Ns^0.5),40];
            p = [1,0.55,1,0.33,2,ceil(0.0625*rndom.Ns^0.5),40];
            
            % Create rndom.filename for each distribtuion
            rndom.filename = sprintf(['D_', char(actual.dist_name),'_T_','%d', '_S_','%d'],i, rndom.Ns);
            
            % file path name
            send_file_name = ['D_',...
                char(actual.dist_name),...
                cpu_type,char(rndom.filename),...
                '.dat'];
            
            % initial details for subsample
            sendFileName1 = ['D_',char(actual.dist_name),cpu_type,char(rndom.filename),'.dat'];
            
            rndom.randomVSactual = "random";
            rndom = dist_list(rndom);
            sample = rndom.rndData;
            % while loop to replace non-finite values with finite values
            
            %                 pd = rndom.distInfo;
            %                 non_finite_vals = sum(~isfinite(sample));
            %                 % check sample
            %                 for inf_index = 1:non_finite_vals
            %                     index = isinf(sample);
            %                     num_inf = non_finite_vals;
            %
            %                     new_data = inf;
            %                     disp(['Old: ', num2str(sample(inf_index))])
            %                     while isinf(new_data)
            %                         new_data = random(pd);
            %                         sample(inf_index) = new_data;
            %                     end
            %                     disp(['New: ', num2str(sample(inf_index))])
            %                 end
            %
            % BOOTSTRAP START -----------------------------------------
            
            %dlmwrite(['sample_',rndom.filename,'.dat'],sample, 'delimiter',' ','precision',12)
            
            % boot strap sample details
            bootrandomNs = floor(percSample*rndom.Ns);
            numSample = numSubs;
            
            subSampPlot = [];
            % loop to get estimate for bootstrapped samples
            for s = 1:numSample
                % Bootstrap
                subSample = datasample(sample,bootrandomNs,'Replace',false);
                send_file_name = ['BS000',num2str(s),'.dat'];
               
                x = sort(subSample);
                
                % track T,BR per trial
                [pList,T,BRlevel,BR0] = block_definition.r_tree(x,actual.min_limit,actual.max_limit,rndom.filename,p);
                
                T_track(k) = T;
                BR0_track(k,i) = BR0;
                BR_track(k,i) = sum(BRlevel{end,1});
                
                
                % BOOTSTRAP END -------------------------------------------
                
            end
            
        end
            rndom.filename = sprintf(['D_', char(actual.dist_name),'_T_','%d'],i);
    end
        
        dist_BR0_track{j} = BR0_track;
        dist_BR_track{j} = BR_track;
       
       
        dlmwrite(['BRTrack_SE_',rndom.filename,'.dat'],horzcat(sample_vec',mean(dist_BR0_track{j},2)), 'delimiter',' ','precision',12)
      
end
    toc
    %profile viewer
    % END OF PROGRAM ==========================================================