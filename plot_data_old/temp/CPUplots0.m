% use when driver.m is a script and not a function
clc;clear all; close all;
% class assignment4
actual = distributions;
actual.generate_data = false;

% SUBSAMPLING PARAMETERS---------------------------------------------------
% percentage of sample used to create subsample
percSample = 1;
% number of subsamples to generate
numSubs = 1;
%--------------------------------------------------------------------------

labels2 = {'$2^{9}$','$2^{10}$','$2^{11}$','$2^{12}$','$2^{13}$','$2^{14}$',...
    '$2^{15}$','$2^{16}$','$2^{17}$','$2^{18}$','$2^{19}$','$2^{20}$'};
labels2 = {'$2^{9}$','$2^{10}$','$2^{11}$','$2^{12}$','$2^{13}$','$2^{14}$',...
    '$2^{15}$','$2^{16}$','$2^{17}$','$2^{18}$'};

legNames = {'Tri-modal Normal', 'Normal(5,1)', 'Uniform(4,8)',...
    'Beta(o.5,1.5)', 'Beta(2,0.5)', 'Beta(0.5,0.5)', ...
    'Generalize-Pareto(2,1,0)', 'Stable(0.5,0.05,1,4)'};

% YexpScale = -3;
YexpScale = -2;

tic
% User Options ============================================================
% script switching board
estimator_call_flag =       true;   %<- true/false call SE on/off
estimator_plot_flag =       true;   %<- true/false plot SE results on/off
data_type_flag =            true;   %<- true/false integer powers of 2/real powers of 2
save_graphics =             true;   %<- true/false save .png of plots on/off
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

publicationQuality();

% Example distribution to test %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
distribution_vector = ["Beta-a0p5-b1p5","Beta-a2-b0p5","Beta-a0p5-b0p5",...
    "Bimodal-Normal","BirnbaumSaunders","Burr",...
    "Exponential","Extreme-Value","Gamma","Generalized-Extreme-Value",...
    "Generalized-Pareto","HalfNormal","Normal","Square-periodic",...
    "tLocationScale","Uniform","Uniform-Mix","Weibull","Chisquare",...
    "InverseGaussian","Trimodal-Normal","Stable",...
    "Stable2","Stable3","Stable1","BirnbaumSaunders-Stable"];

distribution_vector = ["Normal","Uniform","Beta-a0p5-b0p5","Generalized-Pareto","Uniform-Mix","Stable"];
distribution_vector = ["Normal","Trimodal-Normal","Uniform","Beta-a0p5-b0p5","Generalized-Pareto"];
distribution_vector = ["Normal","Generalized-Pareto","Stable"];
distribution_vector = ["Generalized-Pareto","Uniform-Mix","Stable"];
distribution_vector = ["Normal","Stable"];%,"Normal"];
distribution_vector = ["Trimodal-Normal"];%,, "Uniform","Beta-a0p5-b0p5"];%,"Generalized-Pareto","Stable","Beta-a0p5-b0p5","Trimodal-Normal","Uniform"];
distribution_vector = ["Trimodal-Normal","Normal", "Uniform","Beta-a0p5-b1p5","Beta-a2-b0p5","Beta-a0p5-b0p5","Uniform-Mix","Generalized-Pareto","Stable"];
% distribution_vector = ["Trimodal-Normal","Normal", "Uniform","Beta-a0p5-b1p5","Beta-a2-b0p5","Beta-a0p5-b0p5"];
% distribution_vector = ["Trimodal-Normal","Normal", "Uniform","Beta-a0p5-b1p5","Uniform-Mix","Generalized-Pareto","Stable"];
% distribution_vector = ["Generalized-Pareto","Stable"];
% distribution_vector = ["Trimodal-Normal","Normal", "Uniform","Beta-a0p5-b1p5"];
distribution_vector = ["Beta-a2-b0p5","Beta-a0p5-b0p5","Uniform-Mix","Generalized-Pareto","Stable"];
distribution_vector = ["Trimodal-Normal","Normal", "Uniform","Beta-a0p5-b1p5","Beta-a2-b0p5","Beta-a0p5-b0p5","Uniform-Mix","Generalized-Pareto","Stable"];
distribution_vector = ["Trimodal-Normal","Normal", "Uniform","Beta-a0p5-b1p5","Beta-a2-b0p5","Beta-a0p5-b0p5","Generalized-Pareto","Stable"];
distribution_vector = ["Trimodal-Normal","Normal", "Uniform"];
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

fail_nmem = zeros(trials,(max_pow-min_pow)/step,length(distribution_vector));
fail_se = zeros(trials,(max_pow-min_pow)/step,length(distribution_vector));

NMEMtripFlag = 0;

% change code to multithreaded version
for j = 1:length(distribution_vector)

    trialIndex = 0;

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

    sample_data = zeros((max_pow-min_pow)/step,1);
    sample_data_box = [];
    kl_se = [];
    kl_nmem = [];
    js_se = [];
    js_nmem = [];
    MSE_se = [];
    MSE_nmem = [];
    all_se = {};
    all_nmem = {};


    NMEMpdfpass = zeros((max_pow-min_pow + 1)/step,trials);

    for i = 1:trials
        % Create vector of  samples
        sample_vec = misc_functions.sample_pow(min_pow,max_pow,data_type_flag,step);

        actual = actual.dist_list();
        labels = [];

        for k = 1:length(sample_vec)

            NMEMtripFlag = 0;

            rndom.Ns = sample_vec(k);
            realx = linspace(actual.min_limit,actual.max_limit,rndom.Ns);

            % Create rndom.filename for each distribtuion
            NMEMdata_pdf = sprintf(['NMEM_pdf_D_', char(actual.dist_name),'_T_','%d', '_S_','%d','.dat'],i, rndom.Ns);
            NMEMdata_cdf = sprintf(['NMEM_cdf_D_', char(actual.dist_name),'_T_','%d', '_S_','%d','.dat'],i, rndom.Ns);
%             NMEMdata_sqr = sprintf(['NMEM_sqr_D_', char(actual.dist_name),'_T_','%d', '_S_','%d','.dat'],i, rndom.Ns);

            disp([char(actual.dist_name),...
                ', Trial: ',num2str(i),'/', num2str(trials), ...
                ' sample size: ',num2str(sample_vec(k))])

            try
                NMEMpdf = importdata(NMEMdata_pdf);
                NMEMcdf = importdata(NMEMdata_cdf);
%                 NMEMsqr = importdata(NMEMdata_sqr);
                NMEMpdfpdfpass(k) = k;
            catch
                NMEMtripFlag = 1;
            end


            if NMEMtripFlag == 1
                break;
            end

            sample_data_box = [sample_data_box,rndom.Ns];
            sample_data(k) = rndom.Ns;

        end


        if NMEMtripFlag == 1
            %break;
        else
            trialIndex = trialIndex + 1;
        end
    end

    NEMMEcolor = [0, 0.75, 0.75];
    SEColor = [0.8,0,0];

    rndom.dist_list();

%     rndom.filename = ['D_', char(rndom.dist_name),'_T_1'];
    rndom.filename = ['D_', char(rndom.dist_name),'_T_3'];

    try
        fail_rate_nmem = importdata(['Failures_NMEM_',rndom.filename,'.dat']);
    catch
        disp(['Failures_NMEM_',rndom.filename,'.dat'])
        break;
    end
    try
        dist_BR0_track_avg{j} = importdata(['BRTrack_SE_',rndom.filename,'.dat']);
    catch
        disp(['BRTrack_SE_',rndom.filename,'.dat'])
        break;
    end
    try
        cpu_vec_nmem = importdata(['cpu_NMEM_',rndom.filename,'.dat']);
    catch
        disp(['cpu_NMEM_',rndom.filename,'.dat'])
        break;
    end
    try
        fail_rate_se = importdata(['Failures_SE_',rndom.filename,'.dat']);
    catch
        disp(['Failures_SE_',rndom.filename,'.dat'])
        break;
    end
    try
        dist_BR0_track_avg{j} = importdata(['BRTrack_SE_',rndom.filename,'.dat']);
    catch
        disp(['BRTrack_SE_',rndom.filename,'.dat'])
        break;
    end
    try
        cpu_vec_se = importdata(['cpu_SE_',rndom.filename,'.dat']);
    catch
        disp(['cpu_SE_',rndom.filename,'.dat'])
        break;
    end

    cpu.stdTimeNMEMpdf{j} = horzcat(cpu_vec_nmem(:,1),std(cpu_vec_nmem(:,(2:end)),[],2));
    cpu.timeNMEMpdf{j} = horzcat(cpu_vec_nmem(:,1),mean(cpu_vec_nmem(:,(2:end)),2));
    cpu.stdTimeSE{j} = horzcat(cpu_vec_se(:,1),std(cpu_vec_se(:,(2:end)),[],2));
    cpu.timeSE{j} = horzcat(cpu_vec_se(:,1),mean(cpu_vec_se(:,(2:end)),2));
    
    % import failed data
    rndom.filename = ['D_', char(rndom.dist_name),'_T_3'];
    
    
    
    fail_se(:,:,j) = fail_rate_se;
    fail_nmem(:,:,j) = fail_rate_nmem;
end
toc
%
% dist_plot_name = ["(a)","(b)","(c)","(d)","(e)"];
dist_plot_name = distribution_vector;

figure('Name','avg Computation Times: cpu')
hold on
cc=lines(length(distribution_vector));
Markers = {'s','p','d','o','^','x','v','^','>','<'};
% for ii = 1:length(distribution_vector)
%
%     xtemp = [ones(length(cpu.timeNMEMpdf{ii}(:,1)),1),log(cpu.timeNMEMpdf{ii}(:,1))/log(2)];
%     ytemp = log(cpu.timeNMEMpdf{ii}(:,2))/log(2);
%     b = xtemp\ytemp;
%     reg_y = xtemp*b;
%
%     h(ii) = plot(log(cpu.timeNMEMpdf{ii}(:,1))/log(2), log(cpu.timeNMEMpdf{ii}(:,2))/log(2),strcat('-',Markers{ii}),'Color',cc(ii,:),...
%         'MarkerEdgeColor',cc(ii,:),...
%         'MarkerFaceColor',cc(ii,:),...
%         'MarkerSize',5',...
%         'DisplayName',char(dist_plot_name(ii)));
% end
for ii = 1:length(distribution_vector)

    xtemp = [ones(length(cpu.timeSE{ii}(:,1)),1),log(cpu.timeSE{ii}(:,1))/log(2)];
    ytemp = log(cpu.timeSE{ii}(:,2))/log(2);
    b = xtemp\ytemp;
    reg_y = xtemp*b;

    h(ii) = plot(log(cpu.timeSE{ii}(:,1))/log(2), log(cpu.timeSE{ii}(:,2))/log(2),strcat('-',Markers{ii}),'Color',cc(ii,:),...
        'MarkerEdgeColor',cc(ii,:),...
        'MarkerFaceColor',cc(ii,:),...
        'MarkerSize',5',...
        'DisplayName',char(dist_plot_name(ii)));
end
for ii = 1:length(distribution_vector)

    xtemp = [ones(length(cpu.timeNMEMpdf{ii}(:,1)),1),log(cpu.timeNMEMpdf{ii}(:,1))/log(2)];
    ytemp = log(cpu.timeNMEMpdf{ii}(:,2))/log(2);
    b = xtemp\ytemp;
    reg_y = xtemp*b;

    g(ii) = plot(log(cpu.timeNMEMpdf{ii}(:,1))/log(2), log(cpu.timeNMEMpdf{ii}(:,2))/log(2),strcat('-k',Markers{ii}),...
        'MarkerEdgeColor',cc(ii,:),...
        'MarkerFaceColor',cc(ii,:),...
        'MarkerSize',5',...
        'DisplayName',char(dist_plot_name(ii)));
end
ylabel('$log_{2}(time)$ (sec)','Interpreter','latex')
xlabel('$log_{2}(N)$','Interpreter','latex')
% legend([h(1:length(distribution_vector))],'Location','NorthWest')
legend(legNames,'Location','NorthWest')
if save_graphics
    graphics_name = 'avg_cpu_times_NMEMpdf';
    png_file = strcat(graphics_name,'.png');
    saveas(gcf,png_file)
    fig_file = strcat(graphics_name,'.eps');
    print(fig_file,'-depsc')

end

figure('Name','std Computation Times: cpu')
hold on
cc=lines(length(distribution_vector));
Markers = {'s','p','d','o','^','x','v','^','>','<'};
% for ii = 1:length(distribution_vector)
%
%     xtemp = [ones(length(cpu.stdTimeNMEMpdf{ii}(:,1)),1),cpu.stdTimeNMEMpdf{ii}(:,1)];
%     ytemp = cpu.stdTimeNMEMpdf{ii}(:,2);
%     b = xtemp\ytemp;
%     reg_y = xtemp*b;
%
%     h(ii) = plot(cpu.stdTimeNMEMpdf{ii}(:,1), cpu.stdTimeNMEMpdf{ii}(:,2),strcat('-',Markers{ii}),'Color',cc(ii,:),...
%         'MarkerEdgeColor',cc(ii,:),...
%         'MarkerFaceColor',cc(ii,:),...
%         'MarkerSize',5',...
%         'DisplayName',char(dist_plot_name(ii)));
% end
for ii = 1:length(distribution_vector)

    xtemp = [ones(length(cpu.stdTimeSE{ii}(:,1)),1),cpu.stdTimeSE{ii}(:,1)];
    ytemp = cpu.stdTimeSE{ii}(:,2);
    b = xtemp\ytemp;
    reg_y = xtemp*b;

    h(ii) = plot(log(cpu.stdTimeSE{ii}(:,1))/log(2), cpu.stdTimeSE{ii}(:,2),strcat('-',Markers{ii}),'Color',cc(ii,:),...
        'MarkerEdgeColor',cc(ii,:),...
        'MarkerFaceColor',cc(ii,:),...
        'MarkerSize',5',...
        'DisplayName',char(dist_plot_name(ii)));
end
for ii = 1:length(distribution_vector)

    xtemp = [ones(length(cpu.stdTimeNMEMpdf{ii}(:,1)),1),cpu.stdTimeNMEMpdf{ii}(:,1)];
    ytemp = cpu.stdTimeNMEMpdf{ii}(:,2);
    b = xtemp\ytemp;
    reg_y = xtemp*b;

    g(ii) = plot(log(cpu.stdTimeNMEMpdf{ii}(:,1))/log(2), cpu.stdTimeNMEMpdf{ii}(:,2),strcat('-k',Markers{ii}),...
        'MarkerEdgeColor',cc(ii,:),...
        'MarkerFaceColor',cc(ii,:),...
        'MarkerSize',5',...
        'DisplayName',char(dist_plot_name(ii)));
end
ylabel('$\sigma$','Interpreter','latex')
xlabel('$N$','Interpreter','latex')
% legend([h(1:length(distribution_vector))],'Location','NorthWest')
legend(legNames,'Location','NorthWest')
if save_graphics
    graphics_name = 'std_cpu_times_NMEMpdf';
    png_file = strcat(graphics_name,'.png');
    saveas(gcf,png_file)
    fig_file = strcat(graphics_name,'.eps');
    print(fig_file,'-depsc')
end

figure('Name',['Failed pdf'])
cc=lines(length(distribution_vector));
Markers = {'s','p','d','o','^','x','v','^','>','<'};
for j = 1:length(distribution_vector)
    hold on
    plot(fail_nmem(:,1,j),fail_nmem(:,2,j),'-r','DisplayName',char(distribution_vector(j)))
    plot(fail_se(:,1,j),fail_se(:,2,j),'-k','DisplayName',char(distribution_vector(j)))

    h(j) = plot(fail_nmem(:,1,j),fail_nmem(:,2,j),strcat('-k',Markers{j}),...
        'MarkerEdgeColor',cc(j,:),...
        'MarkerFaceColor',cc(j,:),...
        'MarkerSize',5',...
        'DisplayName',char(dist_plot_name(j)));

    g(j) = plot(fail_se(:,1,j),fail_se(:,2,j),strcat('-',Markers{j}),'Color',cc(j,:),...
        'MarkerEdgeColor',cc(j,:),...
        'MarkerFaceColor',cc(j,:),...
        'MarkerSize',5',...
        'DisplayName',char(dist_plot_name(j)));
    
    xlabel('N','Interpreter','latex')
    ylabel('Failure Rate','Interpreter','latex')
    ylim([0,1.125])
    legend([g])
end

%------------------------------------------------------------------------
%{


for j = 1:length(distribution_vector)
    figure('Name',['Failed pdfe_',char(distribution_vector(j))])
    hold on
    plot(sample_vec,sum(fail_nmem(:,:,j),1)/trials,'-r','DisplayName',char(distribution_vector(j)))
    plot(sample_vec,sum(fail_se(:,:,j),1)/trials,'-k','DisplayName',char(distribution_vector(j)))
    xlabel('N','Interpreter','latex')
    ylabel('Failure Rate','Interpreter','latex')
    ylim([0,1.125])
    if save_graphics
        graphics_name = ['failed_Estimates_',...
            char(actual.dist_name),...
            'T_',int2str(trials)];
        png_file = strcat(graphics_name,'.png');
        saveas(gcf,png_file)
        fig_file = strcat(graphics_name,'.eps');
        print(fig_file,'-depsc')
    end
end

% dist_plot_name = ["(a)","(b)","(c)","(d)","(e)"];
dist_plot_name = distribution_vector;

figure('Name','avg Computation Times: cpu')
hold on
cc=lines(length(distribution_vector));
Markers = {'s','p','d','o','^','x','v','^','>','<'};
for ii = 1:length(distribution_vector)

    xtemp = [ones(length(cpu.timeSE{ii}(:,1)),1),log(cpu.timeSE{ii}(:,1))/log(2)];
    ytemp = log(cpu.timeSE{ii}(:,2))/log(2);
    b = xtemp\ytemp;
    reg_y = xtemp*b;

    h(ii) = plot(log(cpu.timeSE{ii}(:,1))/log(2), log(cpu.timeSE{ii}(:,2))/log(2),strcat('-',Markers{ii}),'Color',cc(ii,:),...
        'MarkerEdgeColor',cc(ii,:),...
        'MarkerFaceColor',cc(ii,:),...
        'MarkerSize',5',...
        'DisplayName',char(dist_plot_name(ii)));
end
for ii = 1:length(distribution_vector)

    xtemp = [ones(length(cpu.timeNMEM{ii}(:,1)),1),log(cpu.timeNMEM{ii}(:,1))/log(2)];
    ytemp = log(cpu.timeSE{ii}(:,2))/log(2);
    b = xtemp\ytemp;
    reg_y = xtemp*b;

    g(ii) = plot(log(cpu.timeNMEM{ii}(:,1))/log(2), log(cpu.timeNMEM{ii}(:,2))/log(2),strcat('-k',Markers{ii}),...
        'MarkerEdgeColor',cc(ii,:),...
        'MarkerFaceColor',cc(ii,:),...
        'MarkerSize',5',...
        'DisplayName',char(dist_plot_name(ii)));
end
ylabel('$log_{2}(time)$ (sec)','Interpreter','latex')
xlabel('$log_{2}(Ns)$','Interpreter','latex')
legend([h(1:length(distribution_vector))],'Location','NorthWest')
if save_graphics
    graphics_name = 'avg_cpu_times';
    png_file = strcat(graphics_name,'.png');
    saveas(gcf,png_file)
    fig_file = strcat(graphics_name,'.eps');
    print(fig_file,'-depsc')

end

figure('Name','std Computation Times: cpu')
hold on
cc=lines(length(distribution_vector));
Markers = {'s','p','d','o','^','x','v','^','>','<'};
for ii = 1:length(distribution_vector)

    xtemp = [ones(length(cpu.stdTimeSE{ii}(:,1)),1),cpu.stdTimeSE{ii}(:,1)];
    ytemp = cpu.stdTimeSE{ii}(:,2);
    b = xtemp\ytemp;
    reg_y = xtemp*b;

    h(ii) = plot(cpu.stdTimeSE{ii}(:,1), cpu.stdTimeSE{ii}(:,2),strcat('-',Markers{ii}),'Color',cc(ii,:),...
        'MarkerEdgeColor',cc(ii,:),...
        'MarkerFaceColor',cc(ii,:),...
        'MarkerSize',5',...
        'DisplayName',char(dist_plot_name(ii)));
end
for ii = 1:length(distribution_vector)

    xtemp = [ones(length(cpu.stdTimeNMEM{ii}(:,1)),1),cpu.stdTimeNMEM{ii}(:,2)];
    ytemp = cpu.stdTimeSE{ii}(:,2);
    b = xtemp\ytemp;
    reg_y = xtemp*b;

    g(ii) = plot(cpu.stdTimeNMEM{ii}(:,1), cpu.stdTimeNMEM{ii}(:,2),strcat('-k',Markers{ii}),...
        'MarkerEdgeColor',cc(ii,:),...
        'MarkerFaceColor',cc(ii,:),...
        'MarkerSize',5',...
        'DisplayName',char(dist_plot_name(ii)));
end
ylabel('$log_{2}(time)$ (sec)','Interpreter','latex')
xlabel('$log_{2}(Ns)$','Interpreter','latex')
legend([h(1:length(distribution_vector))],'Location','NorthWest')
if save_graphics
    graphics_name = 'std_cpu_times';
    png_file = strcat(graphics_name,'.png');
    saveas(gcf,png_file)
    fig_file = strcat(graphics_name,'.eps');
    print(fig_file,'-depsc')

end
%}

% END OF PROGRAM ==========================================================
