% function driver(percSample,numSubs)
% ^^ for use with Master_driver.m to different subsampling parameters
profile on

addpath("functions/")
% addpath("compile_nmem/")
addpath("compile_nmem_mv/")

% use when driver.m is a script and not a function
clc;clear; close all;

% error handlingdiary(fullfile('log','error_log_extrema.txt'))
diary on;
filename = ['extrema_script_run-',datestr(datetime(floor(now),'ConvertFrom','datenum')),'.txt'];
full_file = fullfile('log',filename);

fid = fopen(full_file, 'w');
fprintf(fid,['Extrema script started on: ',datestr(datetime(now,'ConvertFrom','datenum')),'/n']);
fclose(fid);

dir_name = fullfile('data','heavy_tails');
status = mkdir(dir_name);

% class assignment
actual = distributions;
actual.generate_data = false;

% YexpScale = -3;
YexpScale = -2;

% User Options ============================================================
% script switching board
estimator_call_flag =       true;   %<- true/false call SE on/off
estimator_plot_flag =       false;   %<- true/false plot SE results on/off
data_type_flag =            true;   %<- true/false integer powers of 2/real powers of 2
save_graphics =             false;   %<- true/false save .png of plots on/off
% rndom data generation parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
max_pow =                   20; %<---- maximum exponent to generate samples
min_pow =                   8; %<---- minimum exponent to generate samples
trials =                    30;  %<--- trials to run to generate heuristics for programs
step =                      1;  %<---- control synthetic rndom samples to skip being created
temp_min_limit =            0; %<---- set upper limit for both
actual.min_limit =          temp_min_limit;  %<--- lower limit to plot
temp_max_limit =            1000; %<---- set upper limit for both
actual.max_limit =          temp_max_limit; %<--- upper limit to plot
x_resolution =              1000;
% changes how data files are accessed when using data generated from the
% Univariant Random Sample Generator available on zmerino's github
cpu_type =                   '\';%<--- '\' or '/' for windows or linux

distribution_vector = ["Beta-a0p5-b1p5","Beta-a2-b0p5","Beta-a0p5-b0p5","Generalized-Pareto","Stable"];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main Function Call Loop used to lable plot figures

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
% fail_nmem = zeros(trials,(max_pow-min_pow)/step,length(distribution_vector));
% fail_nse = zeros(trials,(max_pow-min_pow)/step,length(distribution_vector));
fail_nmem = zeros((max_pow-min_pow+1)/step, trials,length(distribution_vector));
fail_nse = zeros((max_pow-min_pow+1)/step, trials,length(distribution_vector));
lagrange_nse = zeros((max_pow-min_pow+1)/step, trials,length(distribution_vector));
lagrange_nmem = zeros((max_pow-min_pow+1)/step, trials,length(distribution_vector));

global_table = table();
global_tb_nmem = table();
global_tb_nse = table();
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
    sample_vec = utils.sample_pow(min_pow,max_pow,data_type_flag,step);
    cpu_vec_se = zeros(length(sample_vec),trials);
    cpu_vec_nmem = zeros(length(sample_vec),trials);
    kl_vec_se = zeros(length(sample_vec),trials);
    kl_vec_nmem = zeros(length(sample_vec),trials);
    mse_vec_se = zeros(length(sample_vec),trials);
    mse_vec_nmem = zeros(length(sample_vec),trials);
       


        
    for k = 1:length(sample_vec)
        % Initialize cell array to store estimate data
        nse_data = cell(5, trials);
        nmem_data = cell(5, trials);

        for i = 1:trials
        
        actual = actual.dist_list();
            tic  
            
            % initialize to 0. used to count number of failed NMEM pdfe
            fail_flag = 0;
            
            rndom.Ns = sample_vec(k);
            realx = linspace(actual.min_limit,actual.max_limit,rndom.Ns);
            % p-vector definition for Rtree
            p = [1,0.55,1,0.33,2,ceil(0.0625*rndom.Ns^0.5),40];
            
            
            % Create rndom.filename for each distribtuion
            rndom.filename = sprintf(['D_', char(actual.dist_name),'_T_','%d', '_S_','%d'],i, rndom.Ns);
            
            % file path name
            send_file_name = ['D_',...
                char(actual.dist_name),...
                cpu_type,char(rndom.filename),...
                '.dat'];
            
            %==========================================================
            % % % % % % % start of estimates % % % % % % % % %
            %==========================================================
            
            %-- NAP start
            
            sendFileName1 = ['D_',char(actual.dist_name),cpu_type,char(rndom.filename),'.dat'];
            rndom.randomVSactual = "random";
            rndom = dist_list(rndom);
            sample = rndom.rndData;
            
            tintialSE = cputime;

            % nap object instantiation
            nap = NAP;
            nap = nap.stitch(sample);
            
            % extract relevant parameters from object after stich() method
            fail_code = nap.failed;
            SE_x = nap.sx;
            SE_pdf = nap.sPDF;
            SE_cdf = nap.sCDF;
            SE_u = nap.u;
            SE_SQR = nap.sqr;
            nBlocks = nap.nBlocks;
            rndom.Ns = nap.N;
            binrndom.Ns =  nap.binN;
            max_LG = nap.LG_max;
            sum_LG = nap.LG_sum;
            T = nap.T;
            BRlevel = nap.BRlevel;
            BR0 = nap.BR0;

            tcpuSE = cputime-tintialSE;
            
            fail_nse(k,i,j) = fail_code;
                           
            
            %-- NMEM start
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
                        
            %==========================================================
            % % % % % % % % end of estimate % % % % % % % % %
            %==========================================================
            
            % calculate LG multipliers
            lagrange_nse(k,i,j) = sum(max_LG);
            lagrange_nmem(k,i,j) = sum(length(lagrange_multipler));
            
            % calculate kl values
            kl_info_se = actual;
            kl_info_se.x = SE_x;
            kl_info_se = dist_list(kl_info_se);
            
            kl_info_nmem = actual;
            kl_info_nmem.x = x_NMEM';
            kl_info_nmem = dist_list(kl_info_nmem);
            kl_info_nmem;
            %----------------------------------------------------------
            
            SE_test_x = interp1(kl_info_se.x,kl_info_se.pdf_y,SE_x);
            NMEM_test_x = interp1(kl_info_nmem.x,kl_info_nmem.pdf_y,x_NMEM');
            
            if sum(~isfinite(interp1(kl_info_se.x,kl_info_se.pdf_y,SE_x)))...
                    || sum(~isfinite(SE_pdf))... % for some reason indefinite values
                    || sum(~isfinite(pdf_NMEM))...
                    || sum(~isfinite(interp1(kl_info_nmem.x,kl_info_nmem.pdf_y,x_NMEM')))
                disp('non-finite values')
            end
            actual_se_pdf = interp1(kl_info_se.x,kl_info_se.pdf_y,SE_x);
            actual_nmem_pdf = interp1(kl_info_nmem.x,kl_info_nmem.pdf_y,x_NMEM');
            
            % find rand sample to find kl for
            [s_se,idx1] = datasample(SE_pdf,length(x_NMEM));
            [s_nmem,idx2] = datasample(pdf_NMEM',k);
                    
            % read in actual distribution
            rndom.dist_list();
            
            % Create final answer file
            rndom.filename = sprintf(['D_', char(rndom.dist_name),'_T_','%d', '_S_','%d'],i, rndom.Ns);
            
            disp([char(actual.dist_name),...
                ', Trial: ',num2str(i),'/', num2str(trials), ...
                ' sample size: ',num2str(sample_vec(k)), ...
                ' failSE: ', num2str(fail_nse(k,i,j)), ' failNMEM: ', num2str(fail_nmem(k,i,j))])
            
            
            [u_NMEM,sqr_NMEM] = utils.sqr(x_NMEM,pdf_NMEM,sample);


            disp(['Iteration:    j: ', num2str(j),'   k: ', num2str(k), '   i: ', num2str(i)])

            nse_data{1, i} =  SE_x';
            nse_data{2, i} = SE_pdf';
            nse_data{3, i} = SE_cdf';
            nse_data{4, i} = SE_u';
            nse_data{5, i} = SE_SQR';

            nmem_data{1, i} = x_NMEM;
            nmem_data{2, i} = pdf_NMEM;
            nmem_data{3, i} = cdf_NMEM;
            nmem_data{4, i} = u_NMEM;
            nmem_data{5, i} =  sqr_NMEM;

            toc
        end
    
        % Save cell array
        save(fullfile(dir_name,['nse_',num2str(distribution_vector(j)),'_s_',num2str(sample_vec(k)), '.mat']), 'nse_data')
        save(fullfile(dir_name,['nmem_',num2str(distribution_vector(j)),'_s_',num2str(sample_vec(k)),'.mat']), 'nmem_data')

%         fullfile(dir_name,['nse_',num2str(distribution_vector(j)),'_s_',num2str(sample_vec(k)), '.mat'])

    end

end




