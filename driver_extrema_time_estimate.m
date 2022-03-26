 % function driver(percSample,numSubs)
% ^^ for use with Master_driver.m to different subsampling parameters
profile on

% use when driver.m is a script and not a function
clc;clear; close all;

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
max_pow =                   11; %<---- maximum exponent to generate samples
min_pow =                   8; %<---- minimum exponent to generate samples
trials =                    2;  %<--- trials to run to generate heuristics for programs
step =                      1;  %<---- control synthetic rndom samples to skip being created
temp_min_limit =            0; %<---- set upper limit for both
actual.min_limit =          temp_min_limit;  %<--- lower limit to plot
temp_max_limit =            1000; %<---- set upper limit for both
actual.max_limit =          temp_max_limit; %<--- upper limit to plot
x_resolution =              1000;
% changes how data files are accessed when using data generated from the
% Univariant Random Sample Generator available on zmerino's github
cpu_type =                   '\';%<--- '\' or '/' for windows or linux

% empty text file used to track progress
status = mkdir('log');
filename = ['Estimate_Data_Total_Time-', datestr(datetime(floor(now), 'ConvertFrom', 'datenum')), '.txt'];
full_file = fullfile('log', filename);

fid = fopen(full_file, 'w');
fprintf(fid,['Time estimate for estimator data run started on: ', datestr(datetime(now, 'ConvertFrom','datenum')), '\n']);
fprintf(fid,['Sample Range: ', num2str(2^(min_pow)),' - ', num2str(2^(max_pow)), '\n']);
fprintf(fid,['Number of Trials: ', num2str(trials), '\n']);
fclose(fid);

% Example distribution to test %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
distribution_vector = ["Beta-a0p5-b1p5","Beta-a2-b0p5","Beta-a0p5-b0p5",...
    "Bimodal-Normal","BirnbaumSaunders","Burr",...
    "Exponential","Extreme-Value","Gamma","Generalized-Extreme-Value",...
    "Generalized-Pareto","HalfNormal","Normal","Square-periodic",...
    "tLocationScale","Uniform","Uniform-Mix","Weibull","Chisquare",...
    "InverseGaussian","Trimodal-Normal","Stable",...
    "Stable2","Stable3","Stable1","BirnbaumSaunders-Stable"];

distribution_vector = ["Beta-a0p5-b1p5","Beta-a2-b0p5","Beta-a0p5-b0p5","Generalized-Pareto", "Stable"];

distribution = distribution_vector';
names = ["Beta(0.5,1.5)", "Beta(2,0.5)", "Beta(0.5,0.5)", "Generalized-Pareto", "Stable(,,,)"]';
test = table(distribution,names);

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


TotalIterations = length(distribution_vector) * trials * (max_pow - min_pow + 1);

nse_total_time = zeros(1,length(distribution_vector));
nmem_total_time = zeros(1,length(distribution_vector));
total_time = zeros(1,length(distribution_vector));

nse_total_time_per_dist_track = zeros(1,max_pow - min_pow + 1);
nmem_total_time_per_dist_track = zeros(1,max_pow - min_pow + 1);

tTotalStart = tic;
currentIteration = 1;

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
       


        
    for k = 1:length(sample_vec)
        
        nse_x = [];
        nse_pdf = [];
        nse_cdf = [];
        nse_u = [];
        nse_sqr = [];

        nmem_x = [];
        nmem_pdf = [];
        nmem_cdf = [];
        nmem_u = [];
        nmem_sqr = [];

            % Only run one trial to compute estimated run time
        %     for i = 1:trials
                i = 1;
        
        actual = actual.dist_list();
        
            interp_etimate = [];
            estimate_data = {};
            
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
            
            tstart = tic;
            %-- NSE start
            
            sendFileName1 = ['D_',char(actual.dist_name),cpu_type,char(rndom.filename),'.dat'];
            rndom.randomVSactual = "random";
            rndom = dist_list(rndom);
            sample = rndom.rndData;
            
            tintialSE = cputime;
            % mixture model does not return random sample
            [fail_code,SE_x,SE_pdf,SE_cdf,SE_u,SE_SQR,nBlocks,Blacklist,...
                rndom.Ns,binrndom.Ns, max_LG, sum_LG,T,BRlevel,BR0]...
                = stitch_pdf(sample,rndom.filename, actual.min_limit,...
                actual.max_limit,p);

            tcpuSE = cputime-tintialSE;
            
            fail_nse(k,i,j) = fail_code;
            
            tend = toc(tstart);
            PercentComplete = 100 * (currentIteration / TotalIterations);
            currentIteration = currentIteration + 1;

            % Estimate total time per distribution and sample for all
            % trials
            h_conversion = 3600; 
            total_time_per_sample = tend*trials/h_conversion;
            nse_total_time_per_dist_track(1,k) = total_time_per_sample;
            meta_string = ['\n', char(actual.dist_name),...
                ' NSE: Sample size: ', num2str(sample_vec(k)),...
                ' Total Time per Sample: ', num2str(sum(nse_total_time_per_dist_track)), ' hrs'];
                           
            
            disp(meta_string)
        
            fid = fopen(full_file, 'at');
            fprintf(fid, meta_string);
            fclose(fid);


            tstart = tic;               
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
               
            tend = toc(tstart);
            PercentComplete = 100 * (currentIteration / TotalIterations);
            currentIteration = currentIteration + 1;

            % Estimate total time per distribution and sample for all
            % trials
            h_conversion = 3600; 
            total_time_per_sample = tend*trials/h_conversion;
            nmem_total_time_per_dist_track(1,k) = total_time_per_sample;
            meta_string = ['\n', char(actual.dist_name),...
                ' NMEM: Sample size: ', num2str(sample_vec(k)),...
                ' Total Time per Sample: ', num2str(sum(nmem_total_time_per_dist_track)), ' hrs \n'];


            fid = fopen(full_file, 'at');
            fprintf(fid, meta_string);
            fclose(fid);
                        
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
            
            
            [u_NMEM,sqr_NMEM] = misc_functions.sqr(x_NMEM,pdf_NMEM,sample_vec(k));


            nse_x = [nse_x, SE_x'];
            nse_pdf = [nse_pdf, SE_pdf'];
            nse_cdf = [nse_cdf, SE_cdf'];
            nse_u = [nse_u, SE_u'];
            nse_sqr = [nse_sqr, SE_SQR'];

%             nmem_x = [nmem_x, x_NMEM];
%             nmem_pdf = [nmem_pdf, pdf_NMEM];
%             nmem_cdf = [nmem_cdf, cdf_NMEM];
%             nmem_u = [nmem_u, u_NMEM];
%             nmem_sqr = [nmem_sqr, sqr_NMEM];

%         end
    
        writematrix(nse_x,fullfile('data','estimates', ['nse_x_',num2str(distribution_vector(j)),'_s_',num2str(sample_vec(k)),'.txt']),'Delimiter','tab')
        writematrix(nse_pdf,fullfile('data','estimates', ['nse_pdf_',num2str(distribution_vector(j)),'_s_',num2str(sample_vec(k)),'.txt']),'Delimiter','tab')
        writematrix(nse_cdf,fullfile('data','estimates', ['nse_cdf_',num2str(distribution_vector(j)),'_s_',num2str(sample_vec(k)),'.txt']),'Delimiter','tab')
        writematrix(nse_u,fullfile('data','estimates', ['nse_u_',num2str(distribution_vector(j)),'_s_',num2str(sample_vec(k)),'.txt']),'Delimiter','tab')
        writematrix(nse_sqr,fullfile('data','estimates', ['nse_sqr_',num2str(distribution_vector(j)),'_s_',num2str(sample_vec(k)),'.txt']),'Delimiter','tab')

%         writematrix(nmem_x,fullfile('data','estimates', ['nmem_x_',num2str(distribution_vector(j)),'_s_',num2str(sample_vec(k)),'.txt']),'Delimiter','tab')
%         writematrix(nmem_pdf,fullfile('data','estimates', ['nmem_pdf_',num2str(distribution_vector(j)),'_s_',num2str(sample_vec(k)),'.txt']),'Delimiter','tab')
%         writematrix(nmem_cdf,fullfile('data','estimates', ['nmem_cdf_',num2str(distribution_vector(j)),'_s_',num2str(sample_vec(k)),'.txt']),'Delimiter','tab')
%         writematrix(nmem_u,fullfile('data','estimates', ['nmem_u_',num2str(distribution_vector(j)),'_s_',num2str(sample_vec(k)),'.txt']),'Delimiter','tab')
%         writematrix(nmem_sqr,fullfile('data','estimates', ['nmem_sqr_',num2str(distribution_vector(j)),'_s_',num2str(sample_vec(k)),'.txt']),'Delimiter','tab')

    end
        nse_total_time_per_dist = sum(nse_total_time_per_dist_track,2);
        nmem_total_time_per_dist = sum(nmem_total_time_per_dist_track,2);

        nse_total_time(j) = nse_total_time_per_dist;
        nmem_total_time(j) = nmem_total_time_per_dist;
    
        total_time(j) = sum(nse_total_time) + sum(nmem_total_time);
        
        meta_string = ['\n \n',...
            ' NSE: Total Time per Distribution: ', num2str(sum(nse_total_time)), ' hrs',...
            '\n NMEM: Total Time per Distribution: ', num2str(sum(nmem_total_time)), ' hrs'];
        
        disp(meta_string)
    

        fid = fopen(full_file, 'at');
        fprintf(fid, meta_string);
        fprintf(fid,'\n_____________________________________________________');
        fclose(fid);

end
meta_string = ['\nTotal Time of Computation: ', num2str(sum(total_time)), ' hrs'];

disp(meta_string)

fid = fopen(full_file, 'at');
fprintf(fid,'\n');
fprintf(fid,'\n=====================================================');
fprintf(fid, meta_string);
fclose(fid);



