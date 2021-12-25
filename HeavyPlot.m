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
min_pow =                   9; %<---- minimum exponent to generate samples
trials =                    100;  %<--- trials to run to generate heuristics for programs
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
distribution_vector = ["Stable","Generalized-Pareto"];
% distribution_vector = ["Normal"];
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
    %
    %     SEpdfpass = zeros((max_pow-min_pow + 1)/step,trials);
    %     NMEMpdfpass = zeros((max_pow-min_pow + 1)/step,trials);
    SEpdfpass = zeros(trials,(max_pow-min_pow + 1)/step);
    NMEMpdfpass = zeros(trials,(max_pow-min_pow + 1)/step);

    for i = 1:trials
        % Create vector of  samples
        sample_vec = misc_functions.sample_pow(min_pow,max_pow,data_type_flag,step);

        actual = actual.dist_list();
        labels = [];

        %         SEcdf = [];
        %         NMEMcdf = [];

        for k = 1:length(sample_vec)

            rndom.Ns = sample_vec(k);
            realx = linspace(actual.min_limit,actual.max_limit,rndom.Ns);

            % Create rndom.filename for each distribtuion
            SEdata_cdf = sprintf(['SE_cdf_D_', char(actual.dist_name),'_T_','%d', '_S_','%d','.dat'],i, rndom.Ns);
            NMEM_cdf = sprintf(['NMEM_cdf_D_', char(actual.dist_name),'_T_','%d', '_S_','%d','.dat'],i, rndom.Ns);

            disp([char(actual.dist_name),...
                ', Trial: ',num2str(i),'/', num2str(trials), ...
                ' sample size: ',num2str(sample_vec(k))])

            try
                % add if statment
                SEcdf = importdata(SEdata_cdf);
                %                 SEcdf = cell2mat(textscan(fopen(SEdata_cdf),'%f %f'));
                %                 disp('successful Imported')
                %                 disp(['CDF: ',SEdata_cdf])
                %                 disp('test 0')
                %                 size(SEcdf)

                if size(SEcdf,2) == 1
                    SEcdf = [SEcdf,Nan*ones(size(SEcdf,1),1)];
                end

                if any(isnan(SEcdf(:,2))) || size(SEcdf,2) == 1
                    %                     SEtripFlag = 1;
                    %                     break;
                else
                    SEpdfpass(i,k) = 1;
                end
                %                 disp('test 1')
            catch
                %                 disp('Failed Imported')
                %                 disp(['CDF: ',SEdata_cdf])
                SEtripFlag = 1;
                %                 break;
            end

            try
                % add if statment
                NMEMcdf = importdata(NMEM_cdf);
                %                 NMEMcdf = cell2mat(textscan(fopen(NMEM_cdf),'%f %f'));
                %                 disp('successful Imported')
                %                 disp(['CDF: ',NMEM_cdf])
                %                 disp('test 2')

                if size(NMEMcdf,2) == 1
                    NMEMcdf = [NMEMcdf,Nan*ones(size(NMEMcdf,1),1)];
                end

                if any(isnan(NMEMcdf(:,2))) || size(NMEMcdf,2) == 1
                else
                    NMEMpdfpass(i,k) = 1;
                end
                %                 disp('test 3')
            catch
                %                 disp('Failed Imported')
                %                 disp(['CDF: ',NMEM_cdf])
                NMEMtripFlag = 1;
            end
            %             disp('------------')
            % %
            %             NMEMtripFlag
            %             SEtripFlag

            %             if NMEMtripFlag == 1 || SEtripFlag == 1
            %                 break;
            %             end

            labels = [labels, ["2^(" + num2str(min_pow+k-1) + ")"]];

%             SEpdfpass(i,k)
            if SEpdfpass(i,k) == 1
                if any(isnan(SEcdf(:,2))) || size(SEcdf,2) == 1

                else
                    xcdf_SE = SEcdf(:,1);
                    SE_cdf = SEcdf(:,2);

                    all_se{i,k}.xcdf(:,1) = xcdf_SE;
                    all_se{i,k}.cdf(:,1) = SE_cdf;
                end
            end

%             NMEMpdfpass(i,k)
            if NMEMpdfpass(i,k) == 1
                if any(isnan(NMEM_cdf(:,2))) || size(NMEM_cdf,2) == 1
                else
                    xcdf_NMEM = NMEMcdf(:,1);
                    NMEM_cdf = NMEMcdf(:,2);

                    all_nmem{i,k}.xcdf(:,1) = xcdf_NMEM;
                    all_nmem{i,k}.cdf(:,1) = NMEM_cdf;
                end
            end
        end

        if SEtripFlag == 1 || NMEMtripFlag == 1
            %             break;
        else
            trialIndex = trialIndex + 1;
        end

    end

    % SE
    for k = 1:length(sample_vec)
        HTactual = actual;
        figure('Name',['heavy tail compare',char(HTactual.dist_name),...
            'T_',int2str(trials),...
            'S_',int2str(sample_vec(k)),...
            '_',int2str(percSample),...
            '_',int2str(numSubs)])
        hold on
        for i = 1:trials
%             disp('=======')
%             i
%             k
%             SEpdfpass
%             NMEMpdfpass
%             size(all_se,2)
            if SEpdfpass(i,k) ~= 0 %&& size(all_se,2) ~= 1
                basex = 1;%log(exp(1));
                basey = log(1000);
                upperExtend = 1;
                lowerExtend = 1;

                ind = find(all_se{i,k}.xcdf(:,1)>0);
                %                 k
                %                 i
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
                g(1) = plot(xHTse,yHTse,'-b','DisplayName','$\hat{F}(x)_{SE}$');
%                 m = plot(xHTact,yHTact,'-k','DisplayName','$\hat{F}(x)$');
            end
            if NMEMpdfpass(i,k) ~= 0 %&& size(all_nmem,2) ~= 1
                basex = 1;%log(exp(1));
                basey = log(1000);
                upperExtend = 1;
                lowerExtend = 1;

                ind = find(all_nmem{i,k}.xcdf(:,1)>0);

                posXse = all_nmem{i,k}.xcdf(ind,1);
                posCDFse = all_nmem{i,k}.cdf(ind,1);
                xHTse = log(posXse)/basex;
                yHTse = log(1-posCDFse)/basey;

                %             HTactual = actual;

                HTactual.x = posXse;
                HTactual = dist_list(HTactual);
                xHTact = log(HTactual.x)/basex;
                yHTact = log(1-HTactual.cdf_y)/basey;

                %             plot(xHTse,yHTse,'.')
                h(1) = plot(xHTse,yHTse,'-r','DisplayName','$\hat{F}(x)_{NMEM}$');
%                 m = plot(xHTact,yHTact,'-k','DisplayName','$\hat{F}(x)$');
            end
            m(1) = plot(xHTact,yHTact,'-k','DisplayName','$\hat{F}(x)$');
        end
        xlim([lowerExtend*min(xHTact),upperExtend*max(xHTact)]);
        ylim([lowerExtend*min(yHTact),upperExtend*max(yHTact)]);
        ylabel('$\frac{ln(1 - \hat{F}(x))}{3ln(10)}$','Interpreter','latex')
        xlabel('$ln(x)$','Interpreter','latex')
        %         legend('$\hat{f}_{NMEM}(x)$','$\hat{f}_{SE}(x)$','Interpreter','latex','Location','SouthWest')
        legend([h,g,m],'Interpreter','latex','Location','SouthWest')

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
    end

end
toc
% END OF PROGRAM ==========================================================
