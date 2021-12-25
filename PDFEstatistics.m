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

% labels = {'$2^{9}$','$2^{10}$','$2^{11}$','$2^{12}$','$2^{13}$','$2^{14}$','$2^{15}$','$2^{16}$','$2^{17}$'};
% labels = {'$2^{14}$','$2^{15}$','$2^{16}$','$2^{17}$'};
labels = {'$2^{17}$'};
% % labels = {'$2^{9}$','$2^{10}$','$2^{11}$','$2^{12}$','$2^{13}$'};
% labels = {'$2^{10}$','$2^{11}$','$2^{12}$','$2^{13}$','$2^{14}$',...
%     '$2^{15}$','$2^{16}$','$2^{17}$','$2^{18}$','$2^{19}$','$2^{20}$'}; 

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
max_pow =                   17; %<---- maximum exponent to generate samples
min_pow =                   17; %<---- minimum exponent to generate samples
trials =                    1;  %<--- trials to run to generate heuristics for programs
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
distribution_vector = ["Trimodal-Normal","Normal", "Uniform","Beta-a0p5-b1p5","Beta-a2-b0p5","Beta-a0p5-b0p5","Beta-a0p5-b1p5","Beta-a2-b0p5","Uniform-Mix","Generalized-Pareto","Stable"];
% distribution_vector = ["Uniform-Mix","Generalized-Pareto","Stable","Trimodal-Normal","Normal", "Uniform","Beta-a0p5-b1p5","Beta-a2-b0p5","Beta-a0p5-b0p5","Beta-a0p5-b1p5","Beta-a2-b0p5"];
distribution_vector = ["Beta-a0p5-b1p5"];
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

SEtripFlag = 0;
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
    
    for i = 1:trials
        % Create vector of  samples
        sample_vec = misc_functions.sample_pow(min_pow,max_pow,data_type_flag,step);
        
        actual = actual.dist_list();
        
        for k = 1:length(sample_vec)
            
            rndom.Ns = sample_vec(k);
            realx = linspace(actual.min_limit,actual.max_limit,rndom.Ns);
            
            % Create rndom.filename for each distribtuion
            NMEMdata = sprintf(['NMEM_solution_D_', char(actual.dist_name),'_T_','%d', '_S_','%d','.dat'],i, rndom.Ns);
            SEdata = sprintf(['SE_solution_D_', char(actual.dist_name),'_T_','%d', '_S_','%d','.dat'],i, rndom.Ns);
            
            disp([char(actual.dist_name),...
                ', Trial: ',num2str(i),'/', num2str(trials), ...
                ' sample size: ',num2str(sample_vec(k))])%, ...
                %' failSE: ', num2str(fail_se(i,k,j)), ' failNMEM: ', num2str(fail_nmem(i,k,j))])
            
            try
                SE = importdata(SEdata);
            catch
                SEtripFlag = 1;
                break;
            end
            try
                NMEM = importdata(NMEMdata);
            catch
                NMEMtripFlag = 1;
                break;
            end
            
            x_NMEM = NMEM(:,1);
            pdf_NMEM = NMEM(:,2);
            
            x_SE = SE(:,1);
            SE_pdf = SE(:,2);
            
            [x_SE, indexSE] = unique(x_SE);
            [x_NMEM, indexNMEM] = unique(x_NMEM);
            
            % calculate kl values
            kl_info_se = actual;
            kl_info_se.x = x_SE;
            kl_info_se = dist_list(kl_info_se);
            
            kl_info_nmem = actual;
            kl_info_nmem.x = x_NMEM';
            kl_info_nmem = dist_list(kl_info_nmem);
            kl_info_nmem;
            %----------------------------------------------------------
            
            %             [kl_info_se.pdf_y, index] = unique(kl_info_se.pdf_y);
            %
            %             if size(index,1) <= 2
            %                warning('test')
            %             end
            %
            %             [kl_info_se.x, index] = unique(kl_info_se.x);
            %             t4 = size(index)
            
            SE_pdf = SE_pdf(indexSE);
            pdf_NMEM = pdf_NMEM(indexNMEM);
            
            SE_test_x = interp1(kl_info_se.x,kl_info_se.pdf_y,x_SE);
            NMEM_test_x = interp1(kl_info_nmem.x,kl_info_nmem.pdf_y,x_NMEM');
            
            if sum(~isfinite(interp1(kl_info_se.x,kl_info_se.pdf_y,x_SE)))...
                    || sum(~isfinite(SE_pdf))... % for some reason indefinite values
                    || sum(~isfinite(pdf_NMEM))...
                    || sum(~isfinite(interp1(kl_info_nmem.x,kl_info_nmem.pdf_y,x_NMEM')))
                disp('non-finite values')
            end
            actual_se_pdf = interp1(kl_info_se.x,kl_info_se.pdf_y,x_SE);
            actual_nmem_pdf = interp1(kl_info_nmem.x,kl_info_nmem.pdf_y,x_NMEM');
            
            % find rand sample to find kl for
            [s_se,idx1] = datasample(SE_pdf,length(x_NMEM));
            [s_nmem,idx2] = datasample(pdf_NMEM',k);
            
            % kl_se_test = sum(SE_pdf - interp1(actual.x,actual.pdf_y,x))
            % must add code that conditionally accounts for overlaping
            % estimates from either NMEM or SE.
            
            if max(x_NMEM) > max(x_SE) && min(x_SE) > min(x_NMEM)
                kl_se_test = sum(SE_pdf' - interp1(actual.x',actual.pdf_y',x_SE'));
                kl_nmem_test = sum(interp1(x_NMEM,pdf_NMEM,x_SE') - interp1(actual.x',actual.pdf_y',x_SE'));
            elseif max(x_NMEM) < max(x_SE) && min(x_SE) < min(x_NMEM)
                kl_se_test = sum(interp1(x_SE',SE_pdf',x_NMEM) - interp1(actual.x',actual.pdf_y',x_NMEM));
                kl_nmem_test = sum(pdf_NMEM - interp1(actual.x',actual.pdf_y',x_NMEM));
            end
            
            %----------------------------------------------------------
            kl_dist = actual;
            kl_dist.x = x_SE;
            kl_dist = dist_list(kl_dist);
            
            if sum(~isfinite(kl_dist.pdf_y)) ||...
                    sum(~isfinite(SE_pdf)) ||...
                    sum(~isfinite(kl_dist.pdf_y')) ||...
                    sum(~isfinite(pdf_NMEM'))
                warning('the inputs contain non-finite values!')
            end
            
            kl_se = [kl_se,KLDiv(kl_dist.pdf_y', SE_pdf')];
            js_se = [js_se,JSDiv(kl_dist.pdf_y', SE_pdf')];
            
            kl_dist.x = x_NMEM;
            kl_dist = dist_list(kl_dist);
            
            kl_nmem = [kl_nmem,KLDiv(kl_dist.pdf_y', pdf_NMEM')];
            js_nmem = [js_nmem,JSDiv(kl_dist.pdf_y', pdf_NMEM')];
            
            sample_data_box = [sample_data_box,rndom.Ns];
            sample_data(k) = rndom.Ns;
            
            % MSE
            interpSEact = interp1(kl_info_se.x,kl_info_se.pdf_y,x_SE);
            interpNMEMact = interp1(kl_info_nmem.x,kl_info_nmem.pdf_y,x_NMEM');
           
            MSE_se = [MSE_se, sample_vec(k)^(-1)*sum((SE_pdf - interpSEact).^2)];
            MSE_nmem = [MSE_nmem, sample_vec(k)^(-1)*sum((pdf_NMEM' - interpNMEMact).^2)];
            
            % read in actual distribution
            rndom.dist_list();
            
            all_se{i,k}.x(:,1) = x_SE;
            all_se{i,k}.pdf(:,1) = SE_pdf;
            all_nmem{i,k}.x(:,1) = x_NMEM;
            all_nmem{i,k}.pdf(:,1) = pdf_NMEM;
            
        end
        
        %         dlmwrite(['cpu_time_se',rndom.filename,'_',int2str(percSample),'_',int2str(numSubs),'.dat'],cpu_time_se,'delimiter',' ')
        %         dlmwrite(['cpu_time_nmem',rndom.filename,'_',int2str(percSample),'_',int2str(numSubs),'.dat'],cpu_time_nmem,'delimiter',' ')
        
        if SEtripFlag == 1 || NMEMtripFlag == 1
            break;
        else
            trialIndex = trialIndex + 1;
        end
    end
    
    NEMMEcolor = [0, 0.75, 0.75];
    SEColor = [0.8,0,0];
    
    for k = 1:length(sample_vec)
        figure('Name',['Plot all trial estimates: ',num2str(sample_vec(k))])
        hold on
        h(1) = plot(all_se{1,k}.x(:,1), all_se{1,k}.pdf(:,1),'-','Color',SEColor,'DisplayName','$\hat{f}^{(k)}(x)_{SE}$');
        h(2) = plot(all_nmem{1,k}.x(:,1), all_nmem{1,k}.pdf(:,1),'-','Color',NEMMEcolor,'DisplayName','$\hat{f}^{(k)}(x)_{NMEM}$');
        for i = 2:trialIndex
            plot(all_se{i,k}.x(:,1), all_se{i,k}.pdf(:,1),'-','Color',SEColor)
            plot(all_nmem{i,k}.x(:,1), all_nmem{i,k}.pdf(:,1),'-','Color',NEMMEcolor)
            
        end
        h(3) = plot(actual.x,actual.pdf_y,'--k','DisplayName','$f(x)$');
        if max(all_nmem{1,k}.pdf(:,1)) > 2
            ylim([0,6])
        else
            ylim([0,1])
        end
        xlim([actual.min_limit,actual.max_limit])
        legend(h,'Location','NorthWest','Interpreter','latex')
        if save_graphics
            graphics_name = ['pdf_all_',...
                char(actual.dist_name),...
                'S_',num2str(sample_vec(k))];
            png_file = strcat(graphics_name,'.png');
            saveas(gcf,png_file)
            fig_file = strcat(graphics_name,'.fig');
            saveas(gcf,fig_file)
            fig_file = strcat(graphics_name,'.eps');
            print(fig_file,'-depsc')
        end
    end
    
    if max(MSE_se) > max(MSE_nmem)
        yrange = [0 (1.1)*max(MSE_se)];
    else
        yrange = [0 (1.1)*max(MSE_nmem)];
    end
    
    figure('Name','SE MSE for estimates')
    boxplot(MSE_se,sample_data_box,'Labels',labels)
    bp = gca;
    bp.XAxis.TickLabelInterpreter = 'latex';
    xlabel('N','Interpreter','latex')
    ylabel('$L_1$','Interpreter','latex')
    ylim(yrange)
    if save_graphics
        graphics_name = ['MSE_SE_',...
            char(actual.dist_name),...
            'T_',int2str(trials)];
        png_file = strcat(graphics_name,'.png');
        saveas(gcf,png_file)
        fig_file = strcat(graphics_name,'.fig');
        saveas(gcf,fig_file)
        fig_file = strcat(graphics_name,'.eps');
        print(fig_file,'-depsc')
    end
    
    figure('Name','NMEM MSE for estimates')
    boxplot(MSE_nmem,sample_data_box,'Labels',labels)
    bp = gca;
    bp.XAxis.TickLabelInterpreter = 'latex';
    xlabel('N','Interpreter','latex')
    ylabel('$L_1$','Interpreter','latex')
    ylim(yrange)
    if save_graphics
        graphics_name = ['MSE_NMEM_',...
            char(actual.dist_name),...
            'T_',int2str(trials)];
        png_file = strcat(graphics_name,'.png');
        saveas(gcf,png_file)
        fig_file = strcat(graphics_name,'.fig');
        saveas(gcf,fig_file)
        fig_file = strcat(graphics_name,'.eps');
        print(fig_file,'-depsc')
    end
    
    %{
    figure('Name','MSE for estimates')
    set(gcf, 'Position',  [100, 100, 1000, 400])
    subplot(1,2,1)
    boxplot(MSE_se,sample_data_box)
    title('SE MSE')
    xlabel('sample')
    ylabel('MSE')
    ylim(yrange)
    
    subplot(1,2,2)
    boxplot(MSE_nmem,sample_data_box)
    title('NMEM MSE')
    xlabel('sample')
    ylabel('MSE')
    ylim(yrange)
    if save_graphics
        graphics_name = ['MSE_NMEM_',...
            char(actual.dist_name),...
            'T_',int2str(trials)];
        png_file = strcat(graphics_name,'.png');
        saveas(gcf,png_file)
        fig_file = strcat(graphics_name,'.fig');
        saveas(gcf,fig_file)
        fig_file = strcat(graphics_name,'.eps');
        print(fig_file,'-depsc')
    end
    %}
    
    
    if max(kl_se) > max(kl_nmem)
        yrange = [0 (1.1)*max(kl_se)];
    else
        yrange = [0 (1.1)*max(kl_nmem)];
    end
        
    figure('Name','SE KL Box PLot')
    boxplot(kl_se,sample_data_box,'Labels',labels)
    bp = gca;
    bp.XAxis.TickLabelInterpreter = 'latex';
    bp.YAxis.Exponent = YexpScale;
    xlabel('N','Interpreter','latex')
    ylabel('$KL$','Interpreter','latex')
    ylim(yrange)
    if save_graphics
        graphics_name = ['Box_KL_SE_',...
            char(actual.dist_name),...
            'T_',int2str(trials)];
        png_file = strcat(graphics_name,'.png');
        saveas(gcf,png_file)
        fig_file = strcat(graphics_name,'.fig');
        saveas(gcf,fig_file)
        fig_file = strcat(graphics_name,'.eps');
        print(fig_file,'-depsc')
    end
    
    figure('Name','NMEM KL Box PLot')
    boxplot(kl_nmem,sample_data_box,'Labels',labels)
    bp = gca;
    bp.XAxis.TickLabelInterpreter = 'latex';
    bp.YAxis.Exponent = YexpScale;
    xlabel('N','Interpreter','latex')
    ylabel('$KL$','Interpreter','latex')
    ylim(yrange)
    if save_graphics
        graphics_name = ['Box_KL_NMEM_',...
            char(actual.dist_name),...
            'T_',int2str(trials)];
        png_file = strcat(graphics_name,'.png');
        saveas(gcf,png_file)
        fig_file = strcat(graphics_name,'.fig');
        saveas(gcf,fig_file)
        fig_file = strcat(graphics_name,'.eps');
        print(fig_file,'-depsc')
    end
    
    %{
    figure('Name','KL Box PLot')
    set(gcf, 'Position',  [100, 100, 1000, 400])
    subplot(1,2,1)
    boxplot(kl_se,sample_data_box)
    title('SE')
    xlabel('Sample Size')
    ylabel('KL')
    ylim(yrange)
    
    subplot(1,2,2)
    boxplot(kl_nmem,sample_data_box)
    title('NMEM')
    xlabel('Sample Size')
    ylabel('KL')
    ylim(yrange)
    if save_graphics
        graphics_name = ['Box_KL_',...
            char(actual.dist_name),...
            'T_',int2str(trials)];
        png_file = strcat(graphics_name,'.png');
        saveas(gcf,png_file)
        fig_file = strcat(graphics_name,'.fig');
        saveas(gcf,fig_file)
        fig_file = strcat(graphics_name,'.eps');
        print(fig_file,'-depsc')
    end
    %}
    
    
    if max(js_se) > max(js_nmem)
        yrange = [0 (1.1)*max(js_se)];
    else
        yrange = [0 (1.1)*max(js_nmem)];
    end
    
    figure('Name','SE JS Box PLot')
    boxplot(js_se,sample_data_box,'Labels',labels)
    bp = gca;
    bp.XAxis.TickLabelInterpreter = 'latex';
    bp.YAxis.Exponent = YexpScale;
    xlabel('N','Interpreter','latex')
    ylabel('$JS$','Interpreter','latex')
    ylim(yrange)
    if save_graphics
        graphics_name = ['Box_JS_SE_',...
            char(actual.dist_name),...
            'T_',int2str(trials)];
        png_file = strcat(graphics_name,'.png');
        saveas(gcf,png_file)
        fig_file = strcat(graphics_name,'.fig');
        saveas(gcf,fig_file)
        fig_file = strcat(graphics_name,'.eps');
        print(fig_file,'-depsc')
    end
    
    figure('Name','NMEM JS Box PLot')
    boxplot(js_nmem,sample_data_box,'Labels',labels)
    bp = gca;
    bp.XAxis.TickLabelInterpreter = 'latex';
    bp.YAxis.Exponent = YexpScale;
    xlabel('N','Interpreter','latex')
    ylabel('$JS$','Interpreter','latex')
    ylim(yrange)
    if save_graphics
        graphics_name = ['Box_JS_NMEM_',...
            char(actual.dist_name),...
            'T_',int2str(trials)];
        png_file = strcat(graphics_name,'.png');
        saveas(gcf,png_file)
        fig_file = strcat(graphics_name,'.fig');
        saveas(gcf,fig_file)
        fig_file = strcat(graphics_name,'.eps');
        print(fig_file,'-depsc')
    end
   
    %{
    % combine all kl values per sample for a given distribution
    kl_compare = [kl_se, kl_nmem];
    MSE_compare = [MSE_se, MSE_nmem];
    estimator_lables_kl = [repmat('SE  ',[length(kl_se), 1]); repmat('NMEM',[length(kl_nmem), 1])];
    estimator_lables_mse = [repmat('SE  ',[length(kl_se), 1]); repmat('NMEM',[length(kl_nmem), 1])];
    
    figure('Name','KL Compare Box PLot')
    boxplot(kl_compare,estimator_lables_kl,'Labels',labels)
    bp = gca;
    bp.XAxis.TickLabelInterpreter = 'latex';
    bp.YAxis.Exponent = -3;
    xlabel('$N$','Interpreter','latex')
    ylabel('$KL$','Interpreter','latex')
    if save_graphics
        graphics_name = ['KL_Compare_Estimates_',...
            char(actual.dist_name),...
            'T_',int2str(trials)];
        png_file = strcat(graphics_name,'.png');
        saveas(gcf,png_file)
        fig_file = strcat(graphics_name,'.fig');
        saveas(gcf,fig_file)
        fig_file = strcat(graphics_name,'.eps');
        print(fig_file,'-depsc')
    end
    
    figure('Name','MSE Compare Box PLot')
    boxplot(MSE_compare,estimator_lables_mse,'Labels',labels)
    bp = gca;
    bp.XAxis.TickLabelInterpreter = 'latex';
    bp.YAxis.Exponent = -3;
    xlabel('$N$','Interpreter','latex')
    ylabel('$L_2$','Interpreter','latex')
    gca.YAxis.Exponent = 0;
    if save_graphics
        graphics_name = ['MSE_Compare_Estimates_',...
            char(actual.dist_name),...
            'T_',int2str(trials)];
        png_file = strcat(graphics_name,'.png');
        saveas(gcf,png_file)
        fig_file = strcat(graphics_name,'.fig');
        saveas(gcf,fig_file)
        fig_file = strcat(graphics_name,'.eps');
        print(fig_file,'-depsc')
    end
    %}
    
    rndom.dist_list();
    
    % Create final answer file
%     rndom.filename = sprintf(['D_', char(rndom.dist_name),'_T_','%d', '_S_','%d'],i, rndom.Ns);
    rndom.filename = ['D_', char(rndom.dist_name),'_T_1'];
    
    try
        fail_rate_se = importdata(['Failures_SE_',rndom.filename,'.dat']);
    catch
        disp(['Failures_SE_',rndom.filename,'.dat'])
        break;
    end
    try
        fail_rate_nmem = importdata(['Failures_NMEM_',rndom.filename,'.dat']);
    catch
        disp(['Failures_NMEM_',rndom.filename,'.dat'])
        break;
    end
    try
        dist_BR0_track_avg = importdata(['BRTrack_SE_',rndom.filename,'.dat']);
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
    try
        cpu_vec_nmem = importdata(['cpu_NMEM_',rndom.filename,'.dat']);
    catch
        disp(['cpu_NMEM_',rndom.filename,'.dat'])
        break;
    end
    
    size(fail_rate_se)
    size(fail_rate_nmem)
    size(dist_BR0_track_avg)
    size(cpu_vec_se)
    size(cpu_vec_nmem)
    
    cpu.stdTimeSE{j} = horzcat(cpu_vec_se(:,1),std(cpu_vec_se(:,(2:end)),[],2));
    cpu.stdTimeNMEM{j} = horzcat(cpu_vec_nmem(:,1),std(cpu_vec_nmem(:,(2:end)),[],2));
    cpu.timeSE{j} = horzcat(cpu_vec_se(:,1),mean(cpu_vec_se(:,(2:end)),2));
    cpu.timeNMEM{j} = horzcat(cpu_vec_nmem(:,1),mean(cpu_vec_nmem(:,(2:end)),2));
    
end
toc

figure('Name','testBR')
hold on
cc=lines(length(distribution_vector));
for o = 1:length(distribution_vector)
    h(o) = plot(dist_BR0_track_avg(:,1),dist_BR0_track_avg(:,2),'DisplayName',char(distribution_vector(o)),'Color',cc(o,:));
end
%t(1) = plot(sample_vec,T_track,'--k','DisplayName','Threshold');
xlabel('2^x')
ylabel('Threshold')
%legend([h,t])
if save_graphics
    graphics_name = ['BR_track_',...
        char(actual.dist_name),...
        'T_',int2str(trials)];
    png_file = strcat(graphics_name,'.png');
    saveas(gcf,png_file)
    fig_file = strcat(graphics_name,'.eps');
    saveas(gcf,fig_file)
end

for j = 1:length(distribution_vector)
    figure('Name',['Failed pdfe_',char(distribution_vector(j))])
    hold on
    plot(fail_rate_se(:,1),fail_rate_se(:,2),'-r','DisplayName',char(distribution_vector(j)))
    plot(fail_rate_nmem(:,1),fail_rate_nmem(:,2),'-k','DisplayName',char(distribution_vector(j)))
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

% END OF PROGRAM ==========================================================