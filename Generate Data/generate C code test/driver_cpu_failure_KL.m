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

labels2 = {'$2^{9}$','$2^{10}$','$2^{11}$','$2^{12}$','$2^{13}$','$2^{14}$',...
    '$2^{15}$','$2^{16}$','$2^{17}$','$2^{18}$','$2^{19}$','$2^{20}$'};
labels2 = {'$2^{9}$','$2^{10}$','$2^{11}$','$2^{12}$','$2^{13}$','$2^{14}$',...
    '$2^{15}$','$2^{16}$','$2^{17}$','$2^{18}$','$2^{19}$','$2^{20}$'};

% YexpScale = -3;
YexpScale = -2;

% User Options ============================================================
% script switching board
estimator_call_flag =       true;   %<- true/false call SE on/off
estimator_plot_flag =       false;   %<- true/false plot SE results on/off
data_type_flag =            true;   %<- true/false integer powers of 2/real powers of 2
save_graphics =             false;   %<- true/false save .png of plots on/off
% rndom data generation parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
max_pow =                   27; %<---- maximum exponent to generate samples
min_pow =                   19; %<---- minimum exponent to generate samples
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
distribution_vector = ["Normal"];
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
            [fail_code,x,SE_pdf,SE_cdf,SE_u,SE_SQR,nBlocks,Blacklist,rndom.Ns,binrndom.Ns, max_LG, sum_LG,T,BRlevel,BR0]...
                = stitch_pdf(sample,rndom.filename,send_file_name,actual.min_limit,actual.max_limit,p);
            tcpuSE = cputime-tintialSE;
            
            fail_se(i,k,j) = fail_code;
                           
            
            %-- NMEM start
            try
                tintialNMEM = cputime;
                [failed, x_NMEM, pdf_NMEM, cdf_NMEM,sqr_NMEM, ~,lagrange] = EstimatePDF(sample);
                fail_nmem(i,k,j) = 0;
                tcpuNMEM = cputime-tintialNMEM;
            
                n = length(sqr_NMEM);
                dx = 1 / (n + 1);
                u_NMEM = dx:dx:(n * dx);
            catch
                warning(['Problem using function.  Assigning a value of 0.']);
                lagrange = 0;
                x_NMEM = linspace(min(sample),max(sample),length(sample))';
                pdf_NMEM = 0*ones(length(sample),1);
                cdf_NMEM = 0*ones(length(sample),1);
                fail_nmem(i,k,j) = 1;
                % don't record cpu time if estimator failed
                tcpuNMEM = NaN;
            end
                        
            %==========================================================
            % % % % % % % % end of estimate % % % % % % % % %
            %==========================================================
            
            % calculate LG multipliers
            max_LG_se = [max_LG_se, sum(max_LG)];
            max_LG_nmem = [max_LG_nmem, sum(length(lagrange))];
            
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
            
            MSE_se = [MSE_se, sample_vec(k)^(-1)*sum((SE_pdf - interpSEact).^2)];
            
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

            MSE_nmem = [MSE_nmem, sample_vec(k)^(-1)*sum((pdf_NMEM' - interpNMEMpdfact).^2)];

            mse_nmem_test = sample_vec(k)^(-1)*sum((pdf_NMEM' - interpNMEMpdfact).^2);
            
            % end of KL ---------------------------------------------------
            
            
            % read in actual distribution
            rndom.dist_list();
            
            % Create final answer file
            rndom.filename = sprintf(['D_', char(rndom.dist_name),'_T_','%d', '_S_','%d'],i, rndom.Ns);
            
            disp([char(actual.dist_name),...
                ', Trial: ',num2str(i),'/', num2str(trials), ...
                ' sample size: ',num2str(sample_vec(k)), ...
                ' failSE: ', num2str(fail_se(i,k,j)), ' failNMEM: ', num2str(fail_nmem(i,k,j))])
            
            % store time of computation
            cpu_vec_se(k,i) = tcpuSE;
            cpu_vec_nmem(k,i) = tcpuNMEM;
            
            
            sample_data_box = [sample_data_box,rndom.Ns];
            kl_se = [kl_se,kl_se_test];
            kl_nmem = [kl_nmem,kl_nmem_test];
            
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
    
    dlmwrite(['MSE_NSE_',rndom.filename,'.txt'],horzcat(sample_vec',mse_vec_se), 'delimiter',' ','precision',12)
    dlmwrite(['MSE_NMEM_',rndom.filename,'.txt'],horzcat(sample_vec',mse_vec_nmem), 'delimiter',' ','precision',12)
    dlmwrite(['KL_NSE_',rndom.filename,'.txt'],horzcat(sample_vec',kl_vec_se), 'delimiter',' ','precision',12)
    dlmwrite(['KL_NMEM_',rndom.filename,'.txt'],horzcat(sample_vec',kl_vec_nmem), 'delimiter',' ','precision',12)
    dlmwrite(['cpu_SE_',rndom.filename,'.dat'],horzcat(sample_vec',cpu_vec_se), 'delimiter',' ','precision',12)
    dlmwrite(['cpu_NMEM_',rndom.filename,'.dat'],horzcat(sample_vec',cpu_vec_nmem), 'delimiter',' ','precision',12)
    dlmwrite(['Failures_SE_',rndom.filename,'.dat'],horzcat(sample_vec',sum(fail_se(:,:,j),1)'/trials), 'delimiter',' ','precision',12)
    dlmwrite(['Failures_NMEM_',rndom.filename,'.dat'],horzcat(sample_vec',sum(fail_nmem(:,:,j),1)'/trials), 'delimiter',' ','precision',12)


    % test plots
    
%     figure('Name','SE MSE for estimates')
%     boxplot(MSE_se,log(sample_data_box)/log(2))%,'Labels',labels2)
%     bp = gca;
%     bp.XAxis.TickLabelInterpreter = 'latex';
%     xlabel('$log_{2}(N)$','Interpreter','latex')
%     ylabel('$L_1$','Interpreter','latex')
% 
% 
%     figure('Name','SE KL Box PLot')
%     boxplot(kl_se,log(sample_data_box)/log(2)')%,'Labels',labels2)
%     bp = gca;
%     bp.XAxis.TickLabelInterpreter = 'latex';
%     bp.YAxis.Exponent = YexpScale;
%     xlabel('$log_{2}(N)$','Interpreter','latex')
%     ylabel('$KL$','Interpreter','latex')
%      
%     figure('Name','NMEM MSE for estimates')
%     boxplot(MSE_nmem,log(sample_data_box)/log(2))%,'Labels',labels2)
%     bp = gca;
%     bp.XAxis.TickLabelInterpreter = 'latex';
%     xlabel('$log_{2}(N)$','Interpreter','latex')
%     ylabel('$L_1$','Interpreter','latex')
%     
%     figure('Name','NMEM KL Box PLot')
%     boxplot(kl_nmem,log(sample_data_box)/log(2)')%,'Labels',labels2)
%     bp = gca;
%     bp.XAxis.TickLabelInterpreter = 'latex';
%     bp.YAxis.Exponent = YexpScale;
%     xlabel('$log_{2}(N)$','Interpreter','latex')
%     ylabel('$KL$','Interpreter','latex')
    
end
%profile viewer
% END OF PROGRAM ==========================================================