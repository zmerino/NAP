
% use when driver.m is a script and not a function
clc;clear all; close all;

addpath("functions/")
addpath("compile_nmem_mv/")

% class assignment
actual = distributions;
actual.generate_data = false;

% User Options ============================================================
% script switching board
estimator_call_flag =       true;   %<- true/false call SE on/off
estimator_plot_flag =       false;   %<- true/false plot SE results on/off
data_type_flag =            true;   %<- true/false integer powers of 2/real powers of 2
save_graphics =             false;   %<- true/false save .png of plots on/off
% rndom data generation parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
max_pow =                   13; %<---- maximum exponent to generate samples
min_pow =                   13; %<---- minimum exponent to generate samples
trials =                    1   ;  %<--- trials to run to generate heuristics for programs
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
distribution_vector = ["Trimodal-Normal"];

distribution = distribution_vector';
names = ["Normal"]';
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
           
        end
    end


end

figure('Name','NSE PDF')
b = plot(x, SE_pdf);
% b = plot(log(x)/log(2), SE_pdf);
bp = gca;
bp.XAxis.TickLabelInterpreter = 'latex';
xlabel('$log_{2}(N)$','Interpreter','latex')
ylabel('$f(N)$','Interpreter','latex')
legend
saveas(gcf,'nse test.png')

figure('Name','NMEM PDF')
% b = plot(log(x_NMEM)/log(2), pdf_NMEM);
b = plot(x_NMEM, pdf_NMEM);
bp = gca;
bp.XAxis.TickLabelInterpreter = 'latex';
xlabel('$log_{2}(N)$','Interpreter','latex')
ylabel('$f(N)$','Interpreter','latex')
legend
saveas(gcf,'nmem test.png')




