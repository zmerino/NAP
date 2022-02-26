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
    '$2^{15}$','$2^{16}$','$2^{17}$','$2^{18}$','$2^{19}$','$2^{20}$'};

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
distribution_vector = ["Trimodal-Normal","Normal", "Uniform","Beta-a0p5-b1p5","Beta-a2-b0p5","Beta-a0p5-b0p5","Generalized-Pareto","Stable"];
% distribution_vector = ["Trimodal-Normal","Normal", "Uniform","Beta-a0p5-b1p5","Beta-a2-b0p5","Beta-a0p5-b0p5"];
% distribution_vector = ["Trimodal-Normal","Normal", "Uniform","Beta-a0p5-b1p5","Uniform-Mix","Generalized-Pareto","Stable"];
% distribution_vector = ["Generalized-Pareto","Stable"];
% distribution_vector = ["Trimodal-Normal","Normal", "Uniform","Beta-a0p5-b1p5"];
% distribution_vector = ["Stable"];
distribution_vector = ["Trimodal-Normal","Normal","Beta-a0p5-b1p5","Generalized-Pareto"];
distribution_vector = ["Normal"];

legNames = {'Tri-modal Normal', 'Normal(5,1)', 'Uniform(4,8)',...
    'Beta(o.5,1.5)', 'Beta(2,0.5)', 'Beta(0.5,0.5)', ...
    'Generalize-Pareto(2,1,0)', 'Stable(0.5,0.05,1,4)'};

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
                % add if statment
                NMEMpdf = importdata(NMEMdata_pdf);
                disp('successful Imported')
                disp(['PDF: ',NMEMdata_pdf])
                if any(isnan(NMEMdata_pdf(:,2))) || ...
                        sum(~isfinite(NMEMdata_pdf(:,2))) || ...
                            sum(~isfinite(NMEMdata_pdf(:,1)))
                else
                    NMEMpdfpass(k,i) = 1;
                end
            catch
                disp('Failed Imported')
                disp(['PDF: ',NMEMdata_pdf])
                SEtripFlag = 1;
                break;
            end


            labels = [labels, ["2^(" + num2str(min_pow+k-1) + ")"]];

            x_NMEMpdf = NMEMpdf(:,1);
            pdf_NMEMpdf = NMEMpdf(:,2);


            ind = isfinite(pdf_NMEMpdf);

            pdf_NMEMpdf = pdf_NMEMpdf(ind);
            x_NMEMpdf = x_NMEMpdf(ind);

%             [x_NMEMpdf, indexNMEMpdf] = unique(x_NMEMpdf);

            % calculate kl values

            %remove data ouside of (0,1)
            if distribution_vector(j) == "Beta-a0p5-b1p5" || ...
                    distribution_vector(j) == "Beta-a2-b0p5" || ...
                    distribution_vector(j) == "Beta-a0p5-b0p5"
                ind_non_fin = find(1>x_NMEMpdf>0);
                x_NMEMpdf = x_NMEMpdf(ind_non_fin);
                pdf_NMEMpdf = pdf_NMEMpdf(ind_non_fin);
            end


            kl_info_nmem = actual;
            kl_info_nmem.x = x_NMEMpdf';
            kl_info_nmem = dist_list(kl_info_nmem);
            kl_info_nmem;

%             pdf_NMEMpdf = pdf_NMEMpdf(indexNMEMpdf);

            if sum(~isfinite(interp1(kl_info_nmem.x,kl_info_nmem.pdf_y,x_NMEMpdf')))...
                    || sum(~isfinite(pdf_NMEMpdf))
                disp('non-finite values')
            end

            actual_nmem_pdf = interp1(kl_info_nmem.x,kl_info_nmem.pdf_y,x_NMEMpdf);

            %----------------------------------------------------------
            kl_dist = actual;
            kl_dist.x = x_NMEMpdf;

            kl_dist = dist_list(kl_dist);

            if sum(~isfinite(kl_dist.pdf_y')) ||...
                    sum(~isfinite(pdf_NMEMpdf'))
                warning('the inputs contain non-finite values!')
            end

%             kl_nmem = [kl_nmem,KLDiv(kl_dist.pdf_y', pdf_NMEMpdf')];
%             js_nmem = [js_nmem,JSDiv(kl_dist.pdf_y', pdf_NMEMpdf')];
            kl_nmem = [kl_nmem,KLDiv(actual_nmem_pdf', pdf_NMEMpdf')];
%             js_nmem = [js_nmem,JSDiv(actual_nmem_pdf', pdf_NMEMpdf')];


%             kl_dist.x = x_NMEMpdf;
%             kl_dist = dist_list(kl_dist);
            sample_data_box = [sample_data_box,rndom.Ns];
            sample_data(k) = rndom.Ns;

            % MSE

            interpNMEMpdfact = interp1(kl_info_nmem.x,kl_info_nmem.pdf_y,x_NMEMpdf');

            MSE_nmem = [MSE_nmem, sample_vec(k)^(-1)*sum((pdf_NMEMpdf' - interpNMEMpdfact).^2)];

            % read in actual distribution
            rndom.dist_list();

            all_nmem{i,k}.x(:,1) = x_NMEMpdf;
            all_nmem{i,k}.pdf(:,1) = pdf_NMEMpdf;

        end


        if NMEMtripFlag == 1
            %break;
        else
            trialIndex = trialIndex + 1;
        end
    end

    NEMMEcolor = [0, 0.75, 0.75];
    SEColor = [0.8,0,0];


    figure('Name','SE MSE for estimates')
    boxplot(MSE_nmem,log(sample_data_box)/log(2))%,'Labels',labels2)
    bp = gca;
    bp.XAxis.TickLabelInterpreter = 'latex';
    xlabel('$log_{2}(N)$','Interpreter','latex')
    ylabel('$L_1$','Interpreter','latex')
    if save_graphics
        graphics_name = ['MSE_NMEMpdf_',...
            char(actual.dist_name),...
            'T_',int2str(trials)];
        png_file = strcat(graphics_name,'.png');
        saveas(gcf,png_file)
        fig_file = strcat(graphics_name,'.fig');
        saveas(gcf,fig_file)
        fig_file = strcat(graphics_name,'.eps');
        print(fig_file,'-depsc')
    end

    figure('Name','NMEMpdf KL Box PLot')
    boxplot(kl_nmem,log(sample_data_box)/log(2))%,'Labels',labels2)
    bp = gca;
    bp.XAxis.TickLabelInterpreter = 'latex';
    bp.YAxis.Exponent = YexpScale;
    xlabel('$log_{2}(N)$','Interpreter','latex')
    ylabel('$KL$','Interpreter','latex')
    if save_graphics
        graphics_name = ['Box_KL_NMEMpdf_',...
            char(actual.dist_name),...
            'T_',int2str(trials)];
        png_file = strcat(graphics_name,'.png');
        saveas(gcf,png_file)
        fig_file = strcat(graphics_name,'.fig');
        saveas(gcf,fig_file)
        fig_file = strcat(graphics_name,'.eps');
        print(fig_file,'-depsc')
    end

end
toc
% END OF PROGRAM ==========================================================
