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
% labels = {'$2^{17}$'};
% labels = {'$2^{9}$','$2^{10}$','$2^{11}$','$2^{12}$','$2^{13}$'};
labels2 = {'$2^{9}$','$2^{10}$','$2^{11}$'};
labels2 = {'$2^{9}$','$2^{10}$','$2^{11}$','$2^{12}$','$2^{13}$','$2^{14}$',...
    '$2^{15}$','$2^{16}$','$2^{17}$','$2^{18}$','$2^{19}$','$2^{20}$'};

% labels = ['$2^{17}$','$2^{18}$','$2^{19}$','$2^{20}$'];
% labels = ["$2^{17}$","$2^{18}$"];

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
max_pow =                   20; %<---- maximum exponent to generate samples
min_pow =                   11; %<---- minimum exponent to generate samples
trials =                    20;  %<--- trials to run to generate heuristics for programs
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
% distribution_vector = ["Normal", "Uniform"];
% distribution_vector = ["Trimodal-Normal", "Beta-a0p5-b1p5"];
distribution_vector = ["Trimodal-Normal","Normal", "Uniform","Beta-a0p5-b1p5","Beta-a2-b0p5","Beta-a0p5-b0p5","Uniform-Mix","Generalized-Pareto","Stable"];
distribution_vector = ["Trimodal-Normal","Normal", "Uniform","Beta-a0p5-b1p5","Beta-a2-b0p5","Beta-a0p5-b0p5","Generalized-Pareto","Stable"];
distribution_vector = ["Stable"];
% distribution_vector = ["Beta-a0p5-b0p5"];
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
passMat = zeros(trials,max_pow-min_pow+1,length(distribution_vector));

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
    
    SEpdfpass = zeros((max_pow-min_pow + 1)/step,trials);
    
    for i = 1:trials
        % Create vector of  samples
        sample_vec = misc_functions.sample_pow(min_pow,max_pow,data_type_flag,step);
        
        actual = actual.dist_list();
        labels = [];
        
        for k = 1:length(sample_vec)
            
            rndom.Ns = sample_vec(k);
            realx = linspace(actual.min_limit,actual.max_limit,rndom.Ns);
            
            % Create rndom.filename for each distribtuion
            SEdata_pdf = sprintf(['SE_pdf_D_', char(actual.dist_name),'_T_','%d', '_S_','%d','.dat'],i, rndom.Ns);
            SEdata_cdf = sprintf(['SE_cdf_D_', char(actual.dist_name),'_T_','%d', '_S_','%d','.dat'],i, rndom.Ns);
            SEdata_sqr = sprintf(['SE_sqr_D_', char(actual.dist_name),'_T_','%d', '_S_','%d','.dat'],i, rndom.Ns);
%             SEdata_pdf = sprintf(['test_SE_pdf_D_', char(actual.dist_name),'_T_','%d', '_S_','%d','.dat'],i, rndom.Ns);
%             SEdata_cdf = sprintf(['test_SE_cdf_D_', char(actual.dist_name),'_T_','%d', '_S_','%d','.dat'],i, rndom.Ns);
%             SEdata_sqr = sprintf(['test_SE_sqr_D_', char(actual.dist_name),'_T_','%d', '_S_','%d','.dat'],i, rndom.Ns);
            
            disp([char(actual.dist_name),...
                ', Trial: ',num2str(i),'/', num2str(trials), ...
                ' sample size: ',num2str(sample_vec(k))])
            
            try
                SEpdf = importdata(SEdata_pdf);
                SEcdf = importdata(SEdata_cdf);
                SEsqr = importdata(SEdata_sqr);
                
                passMat(i,k,j) = 1;
            catch
                SEtripFlag = 1;
                break;
            end
            
            labels = [labels, ["2^(" + num2str(min_pow+k-1) + ")"]];
            
            xpdf_SE = SEpdf(:,1);
            SE_pdf = SEpdf(:,2);
            xcdf_SE = SEcdf(:,1);
            SE_cdf = SEcdf(:,2);
            xsqr_SE = SEsqr(:,1);
            SE_sqr = SEsqr(:,2);
            
            
            if distribution_vector(j) == "Beta-a0p5-b1p5" || ...
                    distribution_vector(j) == "Beta-a2-b0p5" || ...
                    distribution_vector(j) == "Beta-a0p5-b0p5"
                ind_non_fin = find(1>xpdf_SE>0);
                xpdf_SE = xpdf_SE(ind_non_fin);
                SE_pdf = SE_pdf(ind_non_fin);
            end
            
            [xpdf_SE, indexSE] = unique(xpdf_SE);
            
            % calculate kl values
            kl_info_se = actual;
            kl_info_se.x = xpdf_SE;
            kl_info_se = dist_list(kl_info_se);
            
            %----------------------------------------------------------
            
            SE_pdf = SE_pdf(indexSE);
            
            if sum(~isfinite(interp1(kl_info_se.x,kl_info_se.pdf_y,xpdf_SE)))...
                    || sum(~isfinite(SE_pdf))
                disp('non-finite values')
            end
            
            actual_se_pdf = interp1(kl_info_se.x,kl_info_se.pdf_y,xpdf_SE);
            
            %----------------------------------------------------------
            kl_dist = actual;
            kl_dist.x = xpdf_SE;
            kl_dist = dist_list(kl_dist);
            
            if sum(~isfinite(kl_dist.pdf_y)) ||...
                    sum(~isfinite(SE_pdf))
                warning('the inputs contain non-finite values!')
            end
            
            kl_se = [kl_se,KLDiv(kl_dist.pdf_y', SE_pdf')];
            js_se = [js_se,JSDiv(kl_dist.pdf_y', SE_pdf')];
            
            sample_data_box = [sample_data_box,rndom.Ns];
            sample_data(k) = rndom.Ns;
            
            % MSE
            interpSEact = interp1(kl_info_se.x,kl_info_se.pdf_y,xpdf_SE);
            
            MSE_se = [MSE_se, sample_vec(k)^(-1)*sum((SE_pdf - interpSEact).^2)];
            
            % read in actual distribution
            rndom.dist_list();
            
            all_se{i,k}.x(:,1) = xpdf_SE;
            all_se{i,k}.pdf(:,1) = SE_pdf;
            all_se{i,k}.xcdf(:,1) = xcdf_SE;
            all_se{i,k}.cdf(:,1) = SE_cdf;
            all_se{i,k}.u(:,1) = xsqr_SE;
            all_se{i,k}.sqr(:,1) = SE_sqr;
        end
        
    end
    
    NEMMEcolor = [0, 0.75, 0.75];
    SEColor = [0.8,0,0];
    
    %
    for k = 1:length(length(sum(passMat(:,:,1),1)>0))
        figure('Name',['Plot all trial estimates: ',num2str(sample_vec(k))])
        hold on
        for i = 1:trials
            i
            k
            size(all_se{i,k})
            if ~isempty(all_se{i,k})
                plot(all_se{i,k}.x(:,1), all_se{i,k}.pdf(:,1),'-','DisplayName','$\hat{f}^{(k)}(x)_{SE}$');
            else
                continue;
            end
        end
        h = plot(actual.x,actual.pdf_y,'--k','DisplayName','$f(x)$');
        if max(all_se{1,k}.pdf(:,1)) > 2
            ylim([0,6])
        else
            ylim([0,1])
        end
        xlim([actual.min_limit,actual.max_limit])
        legend(h,'Location','NorthWest','Interpreter','latex')
        if save_graphics
            graphics_name = ['pdf_all_SE_',...
                char(actual.dist_name),...
                'S_',num2str(sample_vec(k))];
            png_file = strcat(graphics_name,'.png');
            saveas(gcf,png_file)
            fig_file = strcat(graphics_name,'.fig');
            saveas(gcf,fig_file)
            fig_file = strcat(graphics_name,'.eps');
            print(fig_file,'-depsc')
        end
        
        %sqr ------------------------------------------------------
        %
        figure('Name',['Plot all sqr estimates: ',...
            int2str(rndom.Ns),...
            ' and ',char(actual.dist_name),...
            ' percSample: ',num2str(percSample),...
            ' numSubs: ',num2str(numSubs)])
        hold on;
        smallN = 256;
        smallN2 = 258;
        graymax = 240;
        range = 0:1/(smallN+1):1;
        muLD = range*(smallN + 1) / (smallN + 1);
        lemonDrop = sqrt(muLD.*(1-muLD)) * 3.4;
        sampleCount2 = (smallN + 2):-1:1;
        colorRange = (255-graymax)*sampleCount2/(smallN + 2);
        base = repmat(graymax, smallN + 2, 1);
        col = (base + colorRange') / 255;
        rgb = [col col col];
        count2 = 1;
        for ii = ceil(smallN2/2):smallN2-1
            ix = [ii ii+1 smallN2-ii smallN2-ii+1];
            fill(range(ix), lemonDrop(ix), rgb(count2, :),'edgecolor','none')
            fill(range(ix), -lemonDrop(ix), rgb(count2, :),'edgecolor','none')
            count2 = count2 + 2;
        end
        plot(muLD,lemonDrop,'k--');
        plot(muLD,-lemonDrop,'k--');
        %         ldg2(1) = plot(u_estimate,sqr_estimate,'-k','DisplayName','BSSE');
        %         ldg2(2) = plot(u_NMEM,sqr_NMEM,'-r','DisplayName','BSSE');
        %         h(1) = plot(all_se{1,k}.u(:,1), all_se{1,k}.sqr(:,1),'-','Color',SEColor,'DisplayName','$\hat{f}^{(k)}(x)_{SE}$');
        %         for i = 2:trialIndex
        %             plot(all_se{i,k}.u(:,1), all_se{i,k}.sqr(:,1),'-','Color',SEColor)
        %         end
        for i = 1:trials
            if ~isempty(all_se{i,k})
                plot(all_se{i,k}.u(:,1), all_se{i,k}.sqr(:,1),'-','DisplayName','$\hat{f}^{(k)}(x)_{SE}$');
            else
                continue;
            end
        end
        xlim([0,1])
        ylabel('$SQR$','Interpreter','latex')
        xlabel('u','Interpreter','latex')
        if save_graphics
            graphics_name = ['SE_SQR_D_',...
                char(actual.dist_name),...
                'T_',int2str(trials),...
                'S_',int2str(sample_vec(k)),...
                '_',int2str(percSample),...
                '_',int2str(numSubs)];
            png_file = strcat(graphics_name,'.png');
            saveas(gcf,png_file)
            fig_file = strcat(graphics_name,'.eps');
            print(fig_file,'-depsc')
        end
        %}
        
        %
        HTactual = actual;
        
        figure('Name',['heavy tail compare',char(HTactual.dist_name),...
            'T_',int2str(trials),...
            'S_',int2str(sample_vec(k)),...
            '_',int2str(percSample),...
            '_',int2str(numSubs)])
        hold on
        for i = 1:length(sum(passMat(:,:,1),1)>0)
            basex = 1;%log(exp(1));
            basey = log(1000);
            upperExtend = 1;
            lowerExtend = 1;
            
            
            if ~isempty(all_se{i,k})
                ind = find(all_se{i,k}.xcdf(:,1)>0);
                
                posXse = all_se{i,k}.xcdf(ind,1);
                posCDFse = all_se{i,k}.cdf(ind,1);
                xHTse = log(posXse)/basex;
                yHTse = log(1-posCDFse)/basey;
                
                %             HTactual = actual;
                
                HTactual.x = posXse;
                HTactual = dist_list(HTactual);
                xHTact = log(HTactual.x)/basex;
                yHTact = log(1-HTactual.cdf_y)/basey;
                
                %             plot(xHTse,yHTse,'.')
                plot(xHTse,yHTse,'-')
                plot(xHTact,yHTact,'-k')
            else
                continue;
            end
        end
        xlim([lowerExtend*min(xHTact),upperExtend*max(xHTact)]);
        ylim([lowerExtend*min(yHTact),upperExtend*max(yHTact)]);
        ylabel('$\frac{ln(1 - \hat{F}(x))}{3ln(10)}$','Interpreter','latex')
        xlabel('$ln(x)$','Interpreter','latex')
        legend('$f(x)$','Interpreter','latex','Location','SouthWest')
        
        graphics_name = ['heavyTail_SE_D_',...
            char(actual.dist_name),...
            'T_',int2str(trials),...
            'S_',int2str(sample_vec(k)),...
            '_',int2str(percSample),...
            '_',int2str(numSubs)];
        png_file = strcat(graphics_name,'.png');
        saveas(gcf,png_file)
        fig_file = strcat(graphics_name,'.eps');
        print(fig_file,'-depsc')
        %}
        
        %{
        basex = 1;%log(exp(1));
            basey = log(1000);
            upperExtend = 1;
            lowerExtend = 1;
            
            posXse = all_se{1,k}.x(:,1)(all_se{1,k}.x(:,1)>=0);
            posCDFse = SE_cdf(x>=0);
            
            xHTse = log(posXse)/basex;
            yHTse = log(1-posCDFse)/basey;
            
            HTactual = actual;
            HTactual.x = posXse;
                HTactual = dist_list(HTactual);
               
            xHTact = log(HTactual.x)/basex;
            
            yHTact = log(1-HTactual.cdf_y)/basey;
            
            figure('Name',['heavy tail compare',char(HTactual.dist_name),...
                'T_',int2str(trials),...
                'S_',int2str(sample_vec(k)),...
                '_',int2str(percSample),...
                '_',int2str(numSubs)])
            hold on
            plot(xHTse,yHTse,'.b')
            plot(xHTact,yHTact,'-k')
            
            xlim([lowerExtend*min(xHTact),upperExtend*max(xHTact)]);
            ylim([lowerExtend*min(yHTact),upperExtend*max(yHTact)]);
            ylabel('$\frac{ln(1 - \hat{F}(x))}{3ln(10)}$','Interpreter','latex')
            xlabel('$ln(x)$','Interpreter','latex')
            legend('$\hat{f}_{SE}(x)$','Interpreter','latex','Location','SouthWest')
            graphics_name = ['heavyTail_D_',...
                char(actual.dist_name),...
                'T_',int2str(trials),...
                'S_',int2str(sample_vec(k)),...
                '_',int2str(percSample),...
                '_',int2str(numSubs)];
            png_file = strcat(graphics_name,'.png');
            saveas(gcf,png_file)
            fig_file = strcat(graphics_name,'.eps');
            print(fig_file,'-depsc')
        %}
    end
    %}
    
    
    figure('Name','SE MSE for estimates')
    boxplot(MSE_se,log(sample_data_box)./log(2))%,'Labels',labels2)
    bp = gca;
    bp.XAxis.TickLabelInterpreter = 'latex';
    xlabel('$log_{2}(N)$','Interpreter','latex')
    ylabel('$L_1$','Interpreter','latex')
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
    
    figure('Name','SE KL Box PLot')
    boxplot(kl_se,log(sample_data_box)./log(2))%,'Labels',labels2)
    bp = gca;
    bp.XAxis.TickLabelInterpreter = 'latex';
    bp.YAxis.Exponent = YexpScale;
    xlabel('$log_{2}(N)$','Interpreter','latex')
    ylabel('$KL$','Interpreter','latex')
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
    %{
%     figure('Name','SE JS Box PLot')
%     boxplot(js_se,sample_data_box,'Labels',labels)
%     bp = gca;
%     bp.XAxis.TickLabelInterpreter = 'latex';
%     bp.YAxis.Exponent = YexpScale;
%     xlabel('N','Interpreter','latex')
%     ylabel('$JS$','Interpreter','latex')
%     ylim(yrange)
%     if save_graphics
%         graphics_name = ['Box_JS_SE_',...
%             char(actual.dist_name),...
%             'T_',int2str(trials)];
%         png_file = strcat(graphics_name,'.png');
%         saveas(gcf,png_file)
%         fig_file = strcat(graphics_name,'.fig');
%         saveas(gcf,fig_file)
%         fig_file = strcat(graphics_name,'.eps');
%         print(fig_file,'-depsc')
%     end
%
%     figure('Name','NMEM JS Box PLot')
%     boxplot(js_nmem,sample_data_box,'Labels',labels)
%     bp = gca;
%     bp.XAxis.TickLabelInterpreter = 'latex';
%     bp.YAxis.Exponent = YexpScale;
%     xlabel('N','Interpreter','latex')
%     ylabel('$JS$','Interpreter','latex')
%     ylim(yrange)
%     if save_graphics
%         graphics_name = ['Box_JS_NMEM_',...
%             char(actual.dist_name),...
%             'T_',int2str(trials)];
%         png_file = strcat(graphics_name,'.png');
%         saveas(gcf,png_file)
%         fig_file = strcat(graphics_name,'.fig');
%         saveas(gcf,fig_file)
%         fig_file = strcat(graphics_name,'.eps');
%         print(fig_file,'-depsc')
%     end
    %}
    
    rndom.dist_list();
    
    % Create final answer file
    %     rndom.filename = sprintf(['D_', char(rndom.dist_name),'_T_','%d', '_S_','%d'],i, rndom.Ns);
    rndom.filename = ['D_', char(rndom.dist_name),'_T_1'];
    rndom.filename = ['D_', char(rndom.dist_name),'_T_20'];
    
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
    
    cpu.stdTimeSE{j} = horzcat(cpu_vec_se(:,1),std(cpu_vec_se(:,(2:end)),[],2));
    cpu.timeSE{j} = horzcat(cpu_vec_se(:,1),mean(cpu_vec_se(:,(2:end)),2));
    
    
    %     figure('Name',['Failed pdfe_',char(distribution_vector(j))])
    %     hold on
    %     plot(sample_vec,fail_rate_se/trials,'-k','DisplayName',char(distribution_vector(j)))
    %     xlabel('N','Interpreter','latex')
    %     ylabel('Failure Rate','Interpreter','latex')
    %     ylim([0,1.125])
    %     if save_graphics
    %         graphics_name = ['failed_Estimates_',...
    %             char(actual.dist_name),...
    %             'T_',int2str(trials)];
    %         png_file = strcat(graphics_name,'.png');
    %         saveas(gcf,png_file)
    %         fig_file = strcat(graphics_name,'.eps');
    %         print(fig_file,'-depsc')
    %     end
    
end
toc
%
figure('Name','testBR')
hold on
cc=lines(length(distribution_vector));
for o = 1:length(distribution_vector)
    h(o) = plot(log(dist_BR0_track_avg{o}(:,1))/log(2),dist_BR0_track_avg{o}(:,2),'DisplayName',char(distribution_vector(o)),'Color',cc(o,:));
end
t(1) = plot(log(sample_vec)/log(2),2.6741.*sample_vec.^(-0.626) + 0.005,'--k','DisplayName','T');
bp = gca;
bp.XAxis.TickLabelInterpreter = 'latex';
% bp.xticklabel = labels;
bp.YAxis.Exponent = YexpScale;
xlabel('$log_{2}(x)$','Interpreter','latex')
% xticklabels(labels)
ylabel('$\xi$','Interpreter','latex')
legend([h])
if save_graphics
    graphics_name = ['BR_track_',...
        char(actual.dist_name),...
        'T_',int2str(trials)];
    png_file = strcat(graphics_name,'.png');
    saveas(gcf,png_file)
    fig_file = strcat(graphics_name,'.eps');
    saveas(gcf,fig_file)
end

% for j = 1:length(distribution_vector)
%     figure('Name',['Failed pdfe_',char(distribution_vector(j))])
%     hold on
%     plot(sample_vec,sum(fail_se(:,:,j),1)/trials,'-k','DisplayName',char(distribution_vector(j)))
%     xlabel('N','Interpreter','latex')
%     ylabel('Failure Rate','Interpreter','latex')
%     ylim([0,1.125])
%     if save_graphics
%         graphics_name = ['failed_Estimates_',...
%             char(actual.dist_name),...
%             'T_',int2str(trials)];
%         png_file = strcat(graphics_name,'.png');
%         saveas(gcf,png_file)
%         fig_file = strcat(graphics_name,'.eps');
%         print(fig_file,'-depsc')
%     end
% end



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
ylabel('$log_{2}(time)$ (sec)','Interpreter','latex')
xlabel('$log_{2}(N)$','Interpreter','latex')
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
ylabel('$\sigma$','Interpreter','latex')
xlabel('$log_{2}(N)$','Interpreter','latex')
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