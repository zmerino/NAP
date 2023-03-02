
% use when driver.m is a script and not a function
clc;clear all; close all;

addpath("functions/")
addpath("cpp_code/")

fid = fopen(full_file, 'w');
fprintf(fid,['Cpu failure distance script started on: ',datestr(datetime(now,'ConvertFrom','datenum')),'/n']);
fclose(fid);

dir_name = fullfile('figures','kde_pdfs');
% dir_name = fullfile('figures','kde_pdfs_mod');
status = mkdir(dir_name);

% class assignment
actual = distributions;
actual.generate_data = false;

% User Options ============================================================
% script switching board
estimator_call_flag =       true;   %<- true/false call SE on/off
estimator_plot_flag =       false;   %<- true/false plot SE results on/off
data_type_flag =            true;   %<- true/false integer powers of 2/real powers of 2
save_figs =                 true;   %<- true/false save .png of plots on/off
% rndom data generation parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
max_pow =                   14; %<---- maximum exponent to generate samples
min_pow =                   14; %<---- minimum exponent to generate samples
trials =                    1   ;  %<--- trials to run to generate heuristics for programs
step =                      4;  %<---- control synthetic rndom samples to skip being created
temp_min_limit =            0; %<---- set upper limit for both
actual.min_limit =          temp_min_limit;  %<--- lower limit to plot
temp_max_limit =            1000; %<---- set upper limit for both
actual.max_limit =          temp_max_limit; %<--- upper limit to plot
x_resolution =              1000;
% changes how data files are accessed when using data generated from the
% Univariant Random Sample Generator available on zmerino's github
cpu_type =                   '\';%<--- '\' or '/' for windows or linux


distribution_vector = ["Trimodal-Normal","Uniform","Normal","Uniform-Mix","Beta-a0p5-b1p5","Beta-a2-b0p5","Beta-a0p5-b0p5","Generalized-Pareto","Stable"];
distribution = distribution_vector';
names = ["Tri-Modal-Normal","Uniform", "Normal","Uniform-Mix", "Beta(0.5,1.5)", "Beta(2,0.5)", "Beta(0.5,0.5)", "Generalized-Pareto", "Stable"];


distribution_vector = ["Trimodal-Normal","Uniform","Normal","Uniform-Mix","Beta-a0p5-b1p5","Beta-a2-b0p5","Beta-a0p5-b0p5", "Generalized-Pareto", "Stable"];
distribution = distribution_vector';
names = ["Tri-Modal-Normal","Uniform", "Normal","Uniform-Mix", "Beta(0.5,1.5)", "Beta(2,0.5)", "Beta(0.5,0.5)", "Generalized-Pareto", "Stable"];


% distribution_vector = ["Normal"];
% distribution = distribution_vector';
% names = ["Normal"];



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

% fail_nmem = zeros((max_pow-min_pow+1)/step, trials,length(distribution_vector));
% fail_nse = zeros((max_pow-min_pow+1)/step, trials,length(distribution_vector));
% lagrange_nse = zeros((max_pow-min_pow+1)/step, trials,length(distribution_vector));
% lagrange_nmem = zeros((max_pow-min_pow+1)/step, trials,length(distribution_vector));
fail_nmem = zeros((max_pow-min_pow)/step, trials,length(distribution_vector));
fail_nse = zeros((max_pow-min_pow)/step, trials,length(distribution_vector));
lagrange_nse = zeros((max_pow-min_pow)/step, trials,length(distribution_vector));
lagrange_nmem = zeros((max_pow-min_pow)/step, trials,length(distribution_vector));

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
    sample_vec = utils.sample_pow(min_pow,max_pow,data_type_flag,step);

    % create empty cells for plotting
    nse_pdfs = cell(length(distribution_vector), trials, length(sample_vec), 2);
    nse_cdfs = cell(length(distribution_vector), trials, length(sample_vec), 2);
    nse_sqrs = cell(length(distribution_vector), trials, length(sample_vec), 2);
    
    nse_kde_pdfs = cell(length(distribution_vector), trials, length(sample_vec), 2);
    nse_kde_cdfs = cell(length(distribution_vector), trials, length(sample_vec), 2);
    nse_kde_sqrs = cell(length(distribution_vector), trials, length(sample_vec), 2);

    kde_pdfs = cell(length(distribution_vector), trials, length(sample_vec), 2);
    kde_cdfs = cell(length(distribution_vector), trials, length(sample_vec), 2);
    kde_sqrs = cell(length(distribution_vector), trials, length(sample_vec), 2);

    nmem_pdfs = cell(length(distribution_vector), trials, length(sample_vec), 2);
    nmem_cdfs = cell(length(distribution_vector), trials, length(sample_vec), 2);
    nmem_sqrs = cell(length(distribution_vector), trials, length(sample_vec), 2);

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

            %-- NAP w/ KDE start

            sendFileName1 = ['D_',char(actual.dist_name),cpu_type,char(rndom.filename),'.dat'];
            rndom.randomVSactual = "random";
            rndom = dist_list(rndom);
            sample = rndom.rndData;

            tintialSE = cputime;

            % nap object instantiation
            nse_kde = NAPkde;
            serial = true;
            nse_kde = nse_kde.stitch(sample, serial);
            
            % extract relevant parameters from object after stich() method
            fail_code = nse_kde.failed;
            x_nse_kde = nse_kde.sx;
            SE_pdf_kde = nse_kde.sPDF;
%             SE_cdf = nse_kde.sCDF;
%             SE_u = nse_kde.u;
%             SE_SQR = nse_kde.sqr;
%             nBlocks = nse_kde.nBlocks;
            rndom.Ns = nse_kde.N;
            binrndom.Ns =  nse_kde.binN;
%             max_LG = nse_kde.LG_max;
%             sum_LG = nse_kde.LG_sum;
%             T = nse_kde.T;
%             BRlevel = nse_kde.BRlevel;
%             BR0 = nse_kde.BR0;
% 
%             tcpuSE = cputime-tintialSE;


            %-- NAP w/o KDE start

            sendFileName1 = ['D_',char(actual.dist_name),cpu_type,char(rndom.filename),'.dat'];
            rndom.randomVSactual = "random";
            rndom = dist_list(rndom);
            sample = rndom.rndData;

            tintialSE = cputime;

            % nap object instantiation
            nap = NAP;
            serial = true;
            nap = nap.stitch(sample, serial);
            
            % extract relevant parameters from object after stich() method
            fail_code = nap.failed;
            x_nse = nap.sx;
            SE_pdf = nap.sPDF;
            SE_cdf = nap.sCDF;
            SE_u = nap.u;
            SE_SQR = nap.sqr;
            nBlocks = nap.nBlocks;
            rndom.Ns = nap.N;
            binrndom.Ns =  nap.binN;
            max_LG = nap.LG_max;
            sum_LG = nap.LG_sum;
            T = nap.T;
            BRlevel = nap.BRlevel;
            BR0 = nap.BR0;

            tcpuSE = cputime-tintialSE;

            disp(['Number of blocks: ', num2str(nBlocks)])

            %-- KDE start
            tintialNMEM = cputime;
            [pdf_KDE,x_KDE] = ksdensity(sample);
            tcpuNMEM = cputime-tintialNMEM;


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
            
            % Stitching ---

%             figure('Name','plt_blockpdf')
%             hold on
%             for b=1:length(nse_kde.block_indx)
%                 plot( nse_kde.blocks_x{nse_kde.block_indx(b)} , nse_kde.blocks_pdf{nse_kde.block_indx(b)}, '-b' )
%             end
%             for b=1:length(nse_kde.block_indx)
%                 plot( nap.blocks_x{nap.block_indx(b)} , nap.blocks_pdf{nap.block_indx(b)}, '-r')
%             end
%             ylabel('$\hat{f}(x)$','Interpreter','latex')
%             xlabel('$x$','Interpreter','latex')
%             if max( nse_kde.blocks_x{nse_kde.block_indx(length(nse_kde.block_indx))}) < 1.1
%                 ylim([0,6])
%             else
%                 ylim([0,1])
%             end

            %==========================================================
            % % % % % % % % end of estimate % % % % % % % % %
            %==========================================================


            % SAVING estimates into table ---------------------------------

            nmem_pdfs{j,i,k,1} = x_NMEM;
            nmem_pdfs{j,i,k,2} = pdf_NMEM;

            nse_pdfs{j,i,k,1} = x_nse;
            nse_pdfs{j,i,k,2} = SE_pdf;

            nse_kde_pdfs{j,i,k,1} = x_nse_kde;
            nse_kde_pdfs{j,i,k,2} = SE_pdf_kde;

            kde_pdfs{j,i,k,1} = x_KDE;
            kde_pdfs{j,i,k,2} = pdf_KDE;

            toc
            disp([char(actual.dist_name),...
                ', Trial: ',num2str(i),'/', num2str(trials), ...
                ' sample size: ',num2str(sample_vec(k))])
        end
    end


    % PDFs ----------------------------------------------------------------
    xl1 = [-0.1, 1.1];
    yl1 = [-0.1, 0.8];
    xl2 = [0, 10];
    yl2 = [-0.1, 6];
    for k = 1:length(sample_vec)

        fig_name = sprintf('pdf_d_%s_%s', convertCharsToStrings(distribution_vector(j)), num2str(sample_vec(k)));
        figure('Name',fig_name)

        subplot(2,2,1)
        hold on;
        for i = 1:trials
            plot(actual.x, actual.pdf_y, '--k', DisplayName='$f(x)$')
            plot(nse_kde_pdfs{j,i,k,1}, nse_kde_pdfs{j,i,k,2}, '-b', DisplayName='$\hat{f}(x)_{NAP \ KDE}$')

            if max(actual.x) > 3
                xlim(xl2)
            else
                xlim(xl1)
            end
            if max(actual.pdf_y) > 1
                ylim(yl2)
            else
                ylim(yl1)
            end

            title('$\hat{f}(x)_{NAP \ KDE}$','interpreter','latex')
            ylabel('f(x)')
        end
        subplot(2,2,2)
        hold on;
        for i = 1:trials
            plot(actual.x, actual.pdf_y, '--k', DisplayName='$f(x)$')
            plot(kde_pdfs{j,i,k,1}, kde_pdfs{j,i,k,2}, '-g', DisplayName='$\hat{f}(x)_{KDE}$')

            if max(actual.x) > 3
                xlim(xl2)
            else
                xlim(xl1)
            end
            if max(actual.pdf_y) > 1
                ylim(yl2)
            else
                ylim(yl1)
            end

            title('$\hat{f}(x)_{KDE}$','interpreter','latex')
        end
        subplot(2,2,3)
        hold on;
        for i = 1:trials
            plot(actual.x, actual.pdf_y, '--k', DisplayName='$f(x)$')
            plot(nse_pdfs{j,i,k,1}, nse_pdfs{j,i,k,2}, '-m', DisplayName='$\hat{f}(x)_{NAP}$')

            if max(actual.x) > 3
                xlim(xl2)
            else
                xlim(xl1)
            end
            if max(actual.pdf_y) > 1
                ylim(yl2)
            else
                ylim(yl1)
            end

            title('$\hat{f}(x)_{NAP}$','interpreter','latex')
            xlabel('x')
            ylabel('f(x)')
        end
        subplot(2,2,4)
        hold on;
        for i = 1:trials
            plot(actual.x, actual.pdf_y, '--k', DisplayName='$f(x)$')
            plot(nmem_pdfs{j,i,k,1}, nmem_pdfs{j,i,k,2}, '-r', DisplayName='$\hat{f}(x)_{NMEM}$')

            if max(actual.x) > 3
                xlim(xl2)
            else
                xlim(xl1)
            end
            if max(actual.pdf_y) > 1
                ylim(yl2)
            else
                ylim(yl1)
            end

            title('$\hat{f}(x)_{NMEM}$','interpreter','latex')
        end
        xlabel('x')
        bp = gca;
        if save_figs
            saveas(bp, fullfile(dir_name, [fig_name, '.png']))
        end
    end

end






