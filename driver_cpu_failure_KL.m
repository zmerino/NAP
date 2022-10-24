
% use when driver.m is a script and not a function
clc;clear all; close all;

addpath("functions/")
% addpath("compile_nmem/")
addpath("compile_nmem_mv/")


% error handling
status = mkdir('log');
diary(fullfile('log','error_log_cpu_failure_distance.txt'))
% empty text file used to track progress
filename = ['cpu_failure_distance_script_run-',datestr(datetime(floor(now),'ConvertFrom','datenum')),'.txt'];
full_file = fullfile('log',filename);

fid = fopen(full_file, 'w');
fprintf(fid,['Cpu failure distance script started on: ',datestr(datetime(now,'ConvertFrom','datenum')),'/n']);
fclose(fid);

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
max_pow =                   22; %<---- maximum exponent to generate samples
min_pow =                   8; %<---- minimum exponent to generate samples
trials =                    500   ;  %<--- trials to run to generate heuristics for programs
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
distribution_vector = ["Generalized-Pareto","Stable","Uniform","Beta-a0p5-b0p5"];
distribution_vector = ["Trimodal-Normal","Normal","Generalized-Pareto","Beta-a0p5-b0p5"];
distribution_vector = ["Trimodal-Normal","Normal","Uniform","Beta-a0p5-b1p5","Beta-a2-b0p5","Beta-a0p5-b0p5","Generalized-Pareto", "Uniform-Mix"];

distribution = distribution_vector';
names = ["Tri-Modal-Normal", "Normal", "Uniform", "Beta(0.5,1.5)", "Beta(2,0.5)", "Beta(0.5,0.5)", "Generalized-Pareto", "Uniform-Mix"]';
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
       
    for i = 1:trials
        
        actual = actual.dist_list();
        
        for k = 1:length(sample_vec)
        
            tic  
            
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
            
            %-- NSE start
            
            sendFileName1 = ['D_',char(actual.dist_name),cpu_type,char(rndom.filename),'.dat'];
            rndom.randomVSactual = "random";
            rndom = dist_list(rndom);
            sample = rndom.rndData;
            
            tintialSE = cputime;
            % mixture model does not return random sample
            [fail_code,x,SE_pdf,SE_cdf,SE_u,SE_SQR,nBlocks,Blacklist,...
                rndom.Ns,binrndom.Ns, max_LG, sum_LG,T,BRlevel,BR0]...
                = stitch_pdf(sample,rndom.filename, actual.min_limit,...
                actual.max_limit,p);

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
            kl_info_se.x = x;
            kl_info_se = dist_list(kl_info_se);
            
            kl_info_nmem = actual;
            kl_info_nmem.x = x_NMEM';
            kl_info_nmem = dist_list(kl_info_nmem);
            kl_info_nmem;
            %----------------------------------------------------------
            
            SE_test_x = interp1(kl_info_se.x,kl_info_se.pdf_y,x);
            NMEM_test_x = interp1(kl_info_nmem.x,kl_info_nmem.pdf_y,x_NMEM');
            
            if sum(~isfinite(interp1(kl_info_se.x,kl_info_se.pdf_y,x)))...
                    || sum(~isfinite(SE_pdf))... % for some reason indefinite values
                    || sum(~isfinite(pdf_NMEM))...
                    || sum(~isfinite(interp1(kl_info_nmem.x,kl_info_nmem.pdf_y,x_NMEM')))
                disp('non-finite values')
            end
            actual_se_pdf = interp1(kl_info_se.x,kl_info_se.pdf_y,x);
            actual_nmem_pdf = interp1(kl_info_nmem.x,kl_info_nmem.pdf_y,x_NMEM');
            
            % find rand sample to find kl for
            [s_se,idx1] = datasample(SE_pdf,length(x_NMEM));
            [s_nmem,idx2] = datasample(pdf_NMEM',k);
                      
            
            % KL calculations ---------------------------------------------
            
            % NSE
            kl_dist = actual;
            kl_dist.x = x;
            kl_dist = dist_list(kl_dist);

            if sum(~isfinite(kl_dist.pdf_y)) ||...
                    sum(~isfinite(SE_pdf))
                warning('NSE: the inputs contain non-finite values!')
            end
                        
            kl_se_test = KLDiv(kl_dist.pdf_y', SE_pdf);
            
            % MSE
            interpSEact = interp1(kl_info_se.x,kl_info_se.pdf_y,x);
                       
            mse_se_test = sample_vec(k)^(-1)*sum((SE_pdf - interpSEact).^2);
            
            % NMEM
            
            kl_dist = actual;
            kl_dist.x = x_NMEM;
            kl_dist = dist_list(kl_dist);

            if sum(~isfinite(kl_dist.pdf_y')) ||...
                    sum(~isfinite(pdf_NMEM'))
                warning('NMEM: the inputs contain non-finite values!')
            end

            kl_info_nmem = actual;
            kl_info_nmem.x = x_NMEM';
            kl_info_nmem = dist_list(kl_info_nmem);
            kl_info_nmem;
            actual_nmem_pdf = interp1(kl_info_nmem.x,kl_info_nmem.pdf_y,x_NMEM);
            
            kl_nmem_test = KLDiv(actual_nmem_pdf', pdf_NMEM');
            
            % MSE
            interpNMEMpdfact = interp1(kl_info_nmem.x,kl_info_nmem.pdf_y,x_NMEM');

            mse_nmem_test = sample_vec(k)^(-1)*sum((pdf_NMEM' - interpNMEMpdfact).^2);
            
            % end of KL ---------------------------------------------------
            
            
            % read in actual distribution
            rndom.dist_list();
            
            % Create final answer file
            rndom.filename = sprintf(['D_', char(rndom.dist_name),'_T_','%d', '_S_','%d'],i, rndom.Ns);
            
            disp([char(actual.dist_name),...
                ', Trial: ',num2str(i),'/', num2str(trials), ...
                ' sample size: ',num2str(sample_vec(k)), ...
                ' failSE: ', num2str(fail_nse(k,i,j)), ' failNMEM: ', num2str(fail_nmem(k,i,j))])
            
            % store time of computation
            cpu_vec_se(k,i) = tcpuSE;
            cpu_vec_nmem(k,i) = tcpuNMEM;
            
            % store KL values
            kl_vec_se(k,i) = kl_se_test;
            kl_vec_nmem(k,i) = kl_nmem_test;
            mse_vec_se(k,i) = mse_se_test;
            mse_vec_nmem(k,i) =  mse_nmem_test;
            
            toc
        end
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

    mse_nse = mse_vec_se;
    mse_nmem = mse_vec_nmem;
    kl_nse = kl_vec_se;
    kl_nmem = kl_vec_nmem;
    cpu_nse = cpu_vec_se;
    cpu_nmem = cpu_vec_nmem;

    temp = vertcat(reshape_groups(sample_vec',mse_vec_se),...
    reshape_groups(sample_vec',mse_vec_nmem));
    sample_power = temp(:,1);
    mse = temp(:,2);

    temp = vertcat(reshape_groups(sample_vec',kl_vec_se),...
    reshape_groups(sample_vec',kl_vec_nmem));
    kl = temp(:,2);
    
    temp = vertcat(reshape_groups(sample_vec',cpu_vec_se),...
    reshape_groups(sample_vec',cpu_vec_nmem));
    cpu_time = temp(:,2);

    distribution = repelem(distribution_vector(j), length(temp(:,2)))';


    name = repelem(names(j), length(temp(:,2)))';

    nse_label = repelem(["NSE"], size(reshape_groups(sample_vec',cpu_vec_se), 1));
    nmem_label = repelem(["NMEM"], size(reshape_groups(sample_vec',cpu_vec_nmem), 1));

    temp = vertcat(reshape_groups(sample_vec',fail_nse(:,:,j)),...
        reshape_groups(sample_vec',fail_nmem(:,:,j)));
    

    fail = temp(:,2);   

    temp = vertcat(reshape_groups(sample_vec',lagrange_nse(:,:,j)),...
        reshape_groups(sample_vec',lagrange_nmem(:,:,j)));
    
    lagrange = temp(:,2);    
    
    estimator = vertcat(nse_label', nmem_label');

    test = table(distribution, name, estimator, sample_power, mse, kl, cpu_time, fail, lagrange);

    global_table = [global_table; test];

end

figure('Name','Advanced Box Plot')
b = boxchart(log(global_table.sample_power)/log(2), global_table.mse, 'GroupByColor',global_table.estimator);
bp = gca;
bp.XAxis.TickLabelInterpreter = 'latex';
xlabel('$log_{2}(N)$','Interpreter','latex')
ylabel('$L_1$','Interpreter','latex')
legend

figure('Name','Advanced Box Plot 2')
b = boxchart(log(global_table.sample_power)/log(2), log(global_table.mse), 'GroupByColor',global_table.distribution);
bp = gca;
bp.XAxis.TickLabelInterpreter = 'latex';
xlabel('$log_{2}(N)$','Interpreter','latex')
ylabel('$L_1$','Interpreter','latex')
legend

disp(global_table)

writetable(global_table,fullfile('data','estimator_meta_data.dat'))



