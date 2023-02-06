% use when driver.m is a script and not a function
clc;clear;close all;

addpath("functions/")

dir_name = fullfile('data','theoretical_threshold');
status = mkdir(dir_name);
dir_name = fullfile('data','estimates');
status = mkdir(dir_name);

% error handling
status = mkdir('log');
diary(fullfile('log','error_log_threoretical_threshold.txt'))
diary on;

% empty text file used to track progress
filename = ['threshold_script_run-',datestr(datetime(floor(now),'ConvertFrom','datenum')),'.txt'];
full_file = fullfile('log',filename);

fid = fopen(full_file, 'w');
fprintf(fid,['Threshold script started on: ',datestr(datetime(now,'ConvertFrom','datenum')),'/n']);
fclose(fid);


% class assignment
actual = distributions;
actual.generate_data = false;
% NOTE: true_quantile is commented out because it is slow with Stable dist

YexpScale = -2;

% User Options ============================================================
data_type_flag =            true;   %<- true/false integer powers of 2/real powers of 2

% rndom data generation parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
max_pow =                   22; %<---- maximum exponent to generate samples
min_pow =                   8; %<---- minimum exponent to generate samples
trials =                    30;  %<--- trials to run to generate heuristics for programs
% rndom data generation parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% max_pow =                   20; %<---- maximum exponent to generate samples
% min_pow =                   8; %<---- minimum exponent to generate samples
% trials =                    3;  %<--- trials to run to generate heuristics for programs

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
%{
 distribution_vector = ["Beta-a0p5-b1p5","Beta-a2-b0p5","Beta-a0p5-b0p5",...
    "Bimodal-Normal","BirnbaumSaunders","Burr",...
    "Exponential","Extreme-Value","Gamma","Generalized-Extreme-Value",...
    "Generalized-Pareto","HalfNormal","Normal","Square-periodic",...
    "tLocationScale","Uniform","Uniform-Mix","Weibull","Chisquare",...
    "InverseGaussian","Trimodal-Normal","Stable",...
    "Stable2","Stable3","Stable1","BirnbaumSaunders-Stable"];
%}
distribution_vector = ["Uniform-Mix","Generalized-Pareto","Stable","Trimodal-Normal","Normal", "Uniform","Beta-a0p5-b1p5","Beta-a2-b0p5","Beta-a0p5-b0p5"];
% distribution_vector = ["Normal","Uniform"];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

% Initialize 3D matirx with zeros that is trials by # of samples sizes by
% # number of distributions.
dist_t_track = zeros(trials,(max_pow-min_pow+1),length(distribution_vector));
dist_br0_track = zeros(trials,(max_pow-min_pow+1),length(distribution_vector));

TotalIterations = length(distribution_vector) * trials * (max_pow - min_pow + 1);

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
    %
    for i = 1:trials
        
        actual = actual.dist_list();
        
        for k = 1:length(sample_vec)
            
            % initialize to 0. used to count number of failed NMEM pdfe
            fail_flag = 0;
            
            rndom.Ns = sample_vec(k);
            realx = linspace(actual.min_limit,actual.max_limit,rndom.Ns);
            % p-vector definition for Rtree
            p = [1,0.55,1,0.33,2,ceil(0.0625*rndom.Ns^0.5),40];
            
            % Create rndom.filename for each distribtuion
            rndom.filename = sprintf(['D_', char(actual.dist_name),...
                '_T_','%d', '_S_','%d'],i, rndom.Ns);
            
            rndom.randomVSactual = "random";
            rndom = dist_list(rndom);
            sample = rndom.rndData;
            
            trig = k+min_pow-1;
%             % skip slow sorts
%             if trig > 19 && i > 100
%                 disp('SKIPPED')
%                 
%                 break
%             end
            % skip for very slow sorts
            if trig > 22 && i > 100
                disp('SKIPPED')
                
                break
            end
            
            tstart = tic;
            
            x = sort(sample);
            
            % Track T,BR per trial
            nse = NSE;
            nse.max_bs = 1e3;
            nse = nse.stitch(sample);

            obt_blocks = blocks;
            obt_blocks = obt_blocks.stitch(sample);
            obt_blocks.sample = x;
            obt_blocks.binNs = nse.binN;
            obt_blocks.Ns = length(x);

            [pList,T,BRlevel,BR0] = obt_blocks.r_tree(obt_blocks);
%             [pList,T,BRlevel,BR0] = block_def.r_tree(x,actual.min_limit,actual.max_limit,rndom.filename,p);
            
            tend = toc(tstart);
            m_conversion = 60;
            h_conversion = 3600;
            
            PercentComplete = 100 * (currentIteration / TotalIterations);
            currentIteration = currentIteration + 1;
            
            % Estiamte total time per distribution and sample for all
            % trials
            total_time = tend*trials;

            if i == 1 && k ==1
                meta_string = ['\n','Estimated time of complettion for ',char(actual.dist_name),...
                    ' distribution of sample size ',num2str(sample_vec(k)),': ',...
                    num2str(total_time/h_conversion),'hrs'];

                disp(meta_string)

                fid = fopen(full_file, 'at');
                fprintf(fid,'\n_____________________________________________________');
                fprintf(fid, meta_string);
                fprintf(fid,'\n=====================================================');
                fclose(fid);

            elseif mod(i,10) ==0 && k == 1
                meta_string = [char(actual.dist_name),...
                    ', Trial: ',num2str(i),'/', num2str(trials), ...
                    ' Sample size: ',num2str(sample_vec(k)),...
                    ' Time: ', num2str(tend), ' sec'];
                
                disp('_____________________________________________________')
                disp(meta_string)
                
                disp('=====================================================')
                disp(['Total Elapsed Time: ', num2str(toc(tTotalStart)), 'sec',...
                    ' Percent Complete: ', num2str(PercentComplete),'%'])
                
                
                fid = fopen(full_file, 'at');
                fprintf(fid,'\n_____________________________________________________');
                fprintf(fid, meta_string);
                fprintf(fid,'\n=====================================================');
                fprintf(fid,['\nTotal Elapsed Time: ', num2str(toc(tTotalStart)),...
                    ', Percent Complete: ', num2str(PercentComplete),'%\n']);
                fclose(fid);
            end
            
            BR0_track(k) = BR0;
            
            disp([char(actual.dist_name),...
                ', Trial: ',num2str(i),'/', num2str(trials), ...
                ' sample size: ',num2str(sample_vec(k))])

        end
        
        dist_t_track(i,:,j) = T_track;
        dist_br0_track(i,:,j) = BR0_track;
        
    end
    
    rndom.filename = sprintf(['D_', char(actual.dist_name),'_T_','%d'],i);
end

% For display purposes

% xT = linspace(2^min_pow,2^max_pow,10000);
% T = 2.6741.*xT.^(-0.626) + 0.005;
% 
% figure('Name','testBR')
% hold on;
% cc=lines(length(distribution_vector));
% for j = 1:length(distribution_vector)
% 
%     xorigin = log(sample_vec)/log(2);
%     yorigin = mean(dist_br0_track(:,:,j),1,'omitnan');
%     stdevorigin = std(dist_br0_track(:,:,j),0,1,'omitnan');
% 
%     x = xorigin;
%     y = yorigin;
%     stdev = stdevorigin;
% 
%     curve1 = y + stdev;
%     curve2 = y - stdev;
% 
%     x2 = [x, fliplr(x)];
%     inBetween = [curve1, fliplr(curve2)];
%     f = fill(x2, inBetween,cc(j,:));
%     set(f,'facealpha',0.25)
%     set(f,'edgealpha',0)
% 
% 
%     h(j) = plot(x,y,'DisplayName',char(distribution_vector(j)),'Color',cc(j,:));
% end
% s = plot(log(xT)/log(2), T, 'k--','DisplayName','\Gamma');
% bp = gca;
% bp.XAxis.TickLabelInterpreter = 'latex';
% % bp.xticklabel = labels;
% bp.YAxis.Exponent = YexpScale;
% xlabel('$log_{2}(x)$','Interpreter','latex')
% ylabel('Threshold','Interpreter','latex')
% legend([s,h])
% 
% figure('Name','testBR 2')
% hold on;
% cc=lines(length(distribution_vector));
% for j = 1:length(distribution_vector)
% 
%     xorigin = log(sample_vec)/log(2);
%     yorigin = mean(dist_br0_track(:,:,j),1,'omitnan');
%     stdevorigin = std(dist_br0_track(:,:,j),0,1,'omitnan');
% 
%     x = xorigin;
%     y = yorigin;
%     stdev = stdevorigin;
% 
%     curve1 = y + stdev;
%     curve2 = y - stdev;
% 
%     x2 = [x, fliplr(x)];
%     inBetween = [curve1, fliplr(curve2)];
%     f = fill(x2, inBetween,cc(j,:));
%     set(f,'facealpha',0.25)
%     set(f,'edgealpha',0)
% 
%     h(j) = plot(x,log(y),'DisplayName',char(distribution_vector(j)),'Color',cc(j,:));
% end
% s = plot(log(xT)/log(2), log(T), 'k--','DisplayName','\Gamma');
% bp = gca;
% bp.XAxis.TickLabelInterpreter = 'latex';
% % bp.xticklabel = labels;
% bp.YAxis.Exponent = YexpScale;
% xlabel('$log_{2}(x)$','Interpreter','latex')
% ylabel('log(Threshold)','Interpreter','latex')
% legend([s,h])

% structure data to be written to CSV file with headers given by the
% variable name

SamplePowers = log(sample_vec')/log(2);

UniformMixMean = mean(dist_br0_track(:,:,1),1,'omitnan')';
UniformMixStdev = std(dist_br0_track(:,:,1),0,1,'omitnan')';

GeneralizedParetoMean =  mean(dist_br0_track(:,:,2),1,'omitnan')';
GeneralizedParetoStdev =  std(dist_br0_track(:,:,2),0,1,'omitnan')';

StableMean = mean(dist_br0_track(:,:,3),1,'omitnan')';
StableStdev = std(dist_br0_track(:,:,3),0,1,'omitnan')';

TrimodalNormalMean = mean(dist_br0_track(:,:,4),1,'omitnan')';
TrimodalNormalStdev = std(dist_br0_track(:,:,4),0,1,'omitnan')';

NormalMean = mean(dist_br0_track(:,:,5),1,'omitnan')';
NormalStdev = std(dist_br0_track(:,:,5),0,1,'omitnan')';

UniformMean = mean(dist_br0_track(:,:,6),1,'omitnan')';
UniformStdev = std(dist_br0_track(:,:,6),0,1,'omitnan')';

BetaA0p5B1p5Mean = mean(dist_br0_track(:,:,7),1,'omitnan')';
BetaA0p5B1p5Stdev = std(dist_br0_track(:,:,7),0,1,'omitnan')';

BetaA2B0p5Mean = mean(dist_br0_track(:,:,8),1,'omitnan')';
BetaA2B0p5Stdev = std(dist_br0_track(:,:,8),0,1,'omitnan')';

BetaA0p5B0p5Mean = mean(dist_br0_track(:,:,9),1,'omitnan')';
BetaA0p5B0p5Stdev = std(dist_br0_track(:,:,9),0,1,'omitnan')';


BR0_data_table = table(SamplePowers, UniformMixMean, UniformMixStdev,...
    GeneralizedParetoMean, GeneralizedParetoStdev, StableMean,...
    StableStdev, TrimodalNormalMean, TrimodalNormalStdev, NormalMean,...
    NormalStdev, UniformMean, UniformStdev, BetaA0p5B1p5Mean, ...
    BetaA0p5B1p5Stdev, BetaA2B0p5Mean, BetaA2B0p5Stdev, ...
    BetaA0p5B0p5Mean, BetaA0p5B0p5Stdev);

writetable(BR0_data_table, fullfile('data','theoretical_threshold','br0_table.dat'))

% END OF PROGRAM ==========================================================
diary on