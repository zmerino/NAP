% use when driver.m is a script and not a function
clc;clear;close all;

% class assignment
actual = distributions;
actual.generate_data = false;
% NOTE: true_quantile is commented out because it is slow with Stable dist

labels2 = {'$2^{9}$','$2^{10}$','$2^{11}$','$2^{12}$','$2^{13}$','$2^{14}$',...
    '$2^{15}$','$2^{16}$','$2^{17}$','$2^{18}$'};

% YexpScale = -3;
YexpScale = -2;

% User Options ============================================================
data_type_flag =            true;   %<- true/false integer powers of 2/real powers of 2

% rndom data generation parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
max_pow =                   12; %<---- maximum exponent to generate samples
min_pow =                   8; %<---- minimum exponent to generate samples
trials =                    10;  %<--- trials to run to generate heuristics for programs
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
% distribution_vector = ["Generalized-Pareto","Stable","Trimodal-Normal","Normal", "Uniform","Beta-a0p5-b1p5","Beta-a2-b0p5","Beta-a0p5-b0p5"];
distribution_vector = ["Generalized-Pareto","Stable","Trimodal-Normal","Normal","Uniform"];
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
    
    % track calculated threshold per distrobution
    T_track = zeros(1,max_pow-min_pow+1);
    BR0_track = zeros(1,max_pow-min_pow+1);
    
    sample_data = zeros(max_pow-min_pow,1);
    sample_data_box = [];
    kl_se = [];
    kl_nmem = [];
    Hellinger_se = [];
    Hellinger_nmem = [];
    max_LG_se = zeros(max_pow-min_pow,trials);
    max_LG_nmem = zeros(max_pow-min_pow,trials);
    cpu_vec_se = zeros(max_pow-min_pow,trials);
    cpu_vec_nmem = zeros(max_pow-min_pow,trials);
    %     dist_BR0_track = cell(max_pow-min_pow);
    %     dist_BR_track = cell(max_pow-min_pow);
    dist_BR0_track = cell(length(distribution_vector));
    dist_BR_track = cell(length(distribution_vector));
    MSE_se = [];
    MSE_nmem = [];
    all_se = cell(max_pow-min_pow,trials);
    all_nmem = cell(max_pow-min_pow,trials);
    cpu_time_se = [];
    cpu_time_nmem = [];
    sample_track = [];
    
    % Create vector of  samples
    sample_vec = misc_functions.sample_pow(min_pow,max_pow,data_type_flag,step);
    BoxCPUtimeSE = zeros((max_pow-min_pow)*trials,2);
    BoxCPUtimeNMEM = zeros((max_pow-min_pow)*trials,2);
    
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
            
            % Create rndom.filename for each distribtuion
            rndom.filename = sprintf(['D_', char(actual.dist_name),'_T_','%d', '_S_','%d'],i, rndom.Ns);
            
            
            % file path name
            send_file_name = ['D_',...
                char(actual.dist_name),...
                cpu_type,char(rndom.filename),...
                '.dat'];
            
            %==========================================================
            % % % % % % % start of estimate % % % % % % % % %
            %==========================================================
            tintialSE = cputime;
            
            % initial details for subsample
            sendFileName1 = ['D_',char(actual.dist_name),cpu_type,char(rndom.filename),'.dat'];
            
            rndom.randomVSactual = "random";
            rndom = dist_list(rndom);
            sample = rndom.rndData;
            
            % track T,BR per trial
            [pList,T,BRlevel,BR0] = block_definition.r_tree(sample,actual.min_limit,actual.max_limit,rndom.filename,p);
            
            BR0_track(k) = BR0;
            T_track(k) = T;
            
            %==========================================================
            % % % % % % % % end of estimate % % % % % % % % %
            %==========================================================
            
        end
        
        dist_t_track(i,:,j) = T_track;
        dist_br0_track(i,:,j) = BR0_track;
        
        disp([char(actual.dist_name),...
            ', Trial: ',num2str(i),'/', num2str(trials), ...
            ' sample size: ',num2str(sample_vec(k))])
    end
    
    rndom.filename = sprintf(['D_', char(actual.dist_name),'_T_','%d'],i);
end



figure('Name','testBR')
hold on;
cc=lines(length(distribution_vector));
plot(log(sample_vec)/log(2), T_track, 'k--','DisplayName','T')
for j = 1:length(distribution_vector)
    for i = 1:trials
        h(j) = plot(log(sample_vec)/log(2),dist_br0_track(i,:,j),'DisplayName',char(distribution_vector(j)),'Color',cc(j,:));
    end
end
bp = gca;
bp.XAxis.TickLabelInterpreter = 'latex';
% bp.xticklabel = labels;
bp.YAxis.Exponent = YexpScale;
xlabel('$log_{2}(x)$','Interpreter','latex')
ylabel('Threshold','Interpreter','latex')
legend(h)
% 
% figure('Name','testBR')
% hold on
% cc=lines(length(distribution_vector));
% for o = 1:length(distribution_vector)
%     h(o) = plot(log(dist_BR0_track_avg{o}(:,1))/log(2),dist_BR0_track_avg{o}(:,2),'DisplayName',char(distribution_vector(o)),'Color',cc(o,:));
% end
% t(1) = plot(log(sample_vec)/log(2),2.6741.*sample_vec.^(-0.626) + 0.005,'--k','DisplayName','T');
% bp = gca;
% bp.XAxis.TickLabelInterpreter = 'latex';
% % bp.xticklabel = labels;
% bp.YAxis.Exponent = YexpScale;
% xlabel('$log_{2}(x)$','Interpreter','latex')
% % xticklabels(labels)
% ylabel('Threshold','Interpreter','latex')
% legend([h])

% dlmwrite(['cpu_SE_',rndom.filename,'.dat'],horzcat(sample_vec',cpu_vec_se), 'delimiter',' ','precision',12)
% dlmwrite(['cpu_NMEM_',rndom.filename,'.dat'],horzcat(sample_vec',cpu_vec_nmem), 'delimiter',' ','precision',12)

% END OF PROGRAM ==========================================================
