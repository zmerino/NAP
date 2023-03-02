
% use when driver.m is a script and not a function
clc;clear all; close all;

addpath("functions/")
addpath("cpp_code/")
% addpath("cpp_code_smooth/")

% figure directory
sub_dir = 'pset7';
fig_dir = fullfile('figures','hyper_parameters',sub_dir);
status = mkdir(fig_dir);

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
max_pow =                   17; %<---- maximum exponent to generate samples
min_pow =                   17; %<---- minimum exponent to generate samples
trials =                    2   ;  %<--- trials to run to generate heuristics for programs
step =                      7;  %<---- control synthetic rndom samples to skip being created
temp_min_limit =            0; %<---- set upper limit for both
actual.min_limit =          temp_min_limit;  %<--- lower limit to plot
temp_max_limit =            1000; %<---- set upper limit for both
actual.max_limit =          temp_max_limit; %<--- upper limit to plot
x_resolution =              1000;
% changes how data files are accessed when using data generated from the
% Univariant Random Sample Generator available on zmerino's github
cpu_type =                   '\';%<--- '\' or '/' for windows or linux


% distribution_vector = ["Trimodal-Normal","Uniform","Normal","Uniform-Mix","Beta-a0p5-b1p5","Beta-a2-b0p5","Beta-a0p5-b0p5","Generalized-Pareto","Stable"];
% distribution = distribution_vector';
% names = ["Tri-Modal-Normal","Uniform", "Normal","Uniform-Mix", "Beta(0.5,1.5)", "Beta(2,0.5)", "Beta(0.5,0.5)", "Generalized-Pareto", "Stable"];

distribution_vector = ["Trimodal-Normal","Uniform","Normal","Uniform-Mix","Beta-a0p5-b1p5","Beta-a2-b0p5","Beta-a0p5-b0p5","Generalized-Pareto","Stable", "Stable2","Stable3"];
distribution = distribution_vector';
names = ["Tri-Modal-Normal","Uniform", "Normal","Uniform-Mix", "Beta(0.5,1.5)", "Beta(2,0.5)", "Beta(0.5,0.5)", "Generalized-Pareto", "Stable", "Stable2","Stable3"];

% 
% distribution_vector = ["Beta-a0p5-b0p5","Normal","Trimodal-Normal","Uniform","Uniform-Mix","Beta-a0p5-b1p5","Beta-a2-b0p5","Generalized-Pareto","Stable", "Stable2","Stable3"];
% distribution = distribution_vector';
% names = ["Beta(0.5,0.5)","Normal","Tri-Modal-Normal","Uniform","Uniform-Mix", "Beta(0.5,1.5)", "Beta(2,0.5)", "Generalized-Pareto", "Stable", "Stable2","Stable3"];


% distribution_vector = ["Beta-a0p5-b1p5","Beta-a2-b0p5","Beta-a0p5-b0p5"];
% distribution = distribution_vector';
% names = ["Beta(0.5,1.5)", "Beta(2,0.5)", "Beta(0.5,0.5)"];

% distribution_vector = ["Normal"];
% distribution = distribution_vector';
% names = ["Normal"];

% distribution_vector = ["Beta-a0p5-b1p5","Beta-a2-b0p5","Beta-a0p5-b0p5",...
%     "Bimodal-Normal","BirnbaumSaunders","Burr",...
%     "Exponential","Extreme-Value","Gamma","Generalized-Extreme-Value",...
%     "Generalized-Pareto","HalfNormal","Normal","Square-periodic",...
%     "tLocationScale","Uniform","Uniform-Mix","Weibull","Chisquare",...
%     "InverseGaussian","Trimodal-Normal","Stable",...
%     "Stable2","Stable3","Stable1","BirnbaumSaunders-Stable"];

% distribution_vector = ["Beta-a0p5-b1p5","Beta-a2-b0p5","Beta-a0p5-b0p5",...
%     "Bimodal-Normal",...
%     "Generalized-Pareto","Normal","Square-periodic",...
%     "Uniform","Uniform-Mix",...
%     "Trimodal-Normal","Stable",...
%     "Stable2","Stable3","Stable1","BirnbaumSaunders-Stable"];
% 
% distribution_vector = ["Beta-a0p5-b1p5","Beta-a2-b0p5","Beta-a0p5-b0p5"];
% distribution_vector = ["Stable", "Stable2","Stable3","Stable1"];


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

    nse_lm = cell(length(sample_vec), trials);
    
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

            %-- NAP start

            sendFileName1 = ['D_',char(actual.dist_name),cpu_type,char(rndom.filename),'.dat'];
            rndom.randomVSactual = "random";
            rndom = dist_list(rndom);
            sample = rndom.rndData;

            tintialSE = cputime;

            % nap object instantiation
            nap = NAP;


            serial = false;

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
            nse_LG = nap.LG;
            nse_LG_val = nap.LG_vals;

            nse_lm{k,i} = nap.LG;

            T = nap.T;
            BRlevel = nap.BRlevel;
            BR0 = nap.BR0;

            figure('Name','plt_blockpdf')
            hold on
            for b=1:length(nap.block_indx)
                plot( nap.blocks_x{nap.block_indx(b)} , nap.blocks_pdf{nap.block_indx(b)} )
            end
            ylabel('$\hat{f}(x)$','Interpreter','latex')
            xlabel('$x$','Interpreter','latex')
            if max( nap.blocks_x{nap.block_indx(length(nap.block_indx))}) < 1.1
                ylim([0,6])
            else
                ylim([0,1])
            end

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


            % SAVING estimates into table ---------------------------------

            nmem_pdfs{j,i,k,1} = x_NMEM;
            nmem_pdfs{j,i,k,2} = pdf_NMEM;

            nse_pdfs{j,i,k,1} = x_nse;
            nse_pdfs{j,i,k,2} = SE_pdf;


            nmem_cdfs{j,i,k,1} = x_NMEM;
            nmem_cdfs{j,i,k,2} = cdf_NMEM;

            nse_cdfs{j,i,k,1} = x_nse;
            nse_cdfs{j,i,k,2} = SE_cdf;


            nmem_sqrs{j,i,k,1} = u_NMEM;
            nmem_sqrs{j,i,k,2} = sqr_NMEM;

            nse_sqrs{j,i,k,1} = SE_u;
            nse_sqrs{j,i,k,2} = SE_SQR;

            % -------------------------------------------------------------

%             figure()
%             plot(x_NMEM,pdf_NMEM, DisplayName='pdf nmem')
%             xlabel('x')
%             ylabel('f(x)')
% 
%             figure()
%             plot(x_nse,SE_pdf, DisplayName='pdf nap')
%             xlabel('x')
%             ylabel('f(x)')
% 
%             figure()
%             plot(x_NMEM,cdf_NMEM, DisplayName='cdf nmem')
%             xlabel('x')
%             ylabel('F(x)')
% 
%             figure()
%             plot(x_nse,SE_cdf, DisplayName='cdf nap')
%             xlabel('x')
%             ylabel('F(x)')
% 
% 
%             figure('Name','SQR')
%             hold on;
%             smallN = 256;
%             smallN2 = 258;
%             graymax = 240;
%             range = 0:1/(smallN+1):1;
%             muLD = range*(smallN + 1) / (smallN + 1);
%             lemonDrop = sqrt(muLD.*(1-muLD)) * 3.4;
%             sampleCount2 = (smallN + 2):-1:1;
%             colorRange = (255-graymax)*sampleCount2/(smallN + 2);
%             base = repmat(graymax, smallN + 2, 1);
%             col = (base + colorRange') / 255;
%             rgb = [col col col];
%             count2 = 1;
%             for ii = ceil(smallN2/2):smallN2-1
%                 ix = [ii ii+1 smallN2-ii smallN2-ii+1];
%                 fill(range(ix), lemonDrop(ix), rgb(count2, :),'edgecolor','none')
%                 fill(range(ix), -lemonDrop(ix), rgb(count2, :),'edgecolor','none')
%                 count2 = count2 + 2;
%             end
%             plot(muLD,lemonDrop,'k--');
%             plot(muLD,-lemonDrop,'k--');
%             %             h = plot(nap(d_indx,n_indx).u, nap(d_indx,n_indx).sqr, 'DisplayName', '$\hat{sqr}_i(x)$');
%             nse_h = plot(SE_u, SE_SQR,'-r', 'DisplayName', '$\hat{sqr}_i^{NAP}(x)$');
%             nmem_h = plot(u_NMEM, sqr_NMEM,'-b', 'DisplayName', '$\hat{sqr}_i^{NMEM}(x)$');
%             bp = gca;
%             xlabel('$x$','Interpreter','latex')
%             ylabel('$\hat{SQR}(x)$','Interpreter','latex')
%             %             legend([nse_h(1),nmem_h(1), g],'Interpreter','latex')
%             xlabel('$x$','Interpreter','latex')
%             ylabel('$sqr(x)$','Interpreter','latex')

            toc
            disp([char(actual.dist_name),...
                ', Trial: ',num2str(i),'/', num2str(trials), ...
                ' sample size: ',num2str(sample_vec(k)), ...
                ' failSE: ', num2str(fail_nse(k,i,j)), ' failNMEM: ', ...
                num2str(fail_nmem(k,i,j))])
        end
    end

    for k = 1:length(sample_vec)
        fig_name = [sub_dir,'LG_per_block_pdf_d_',convertStringsToChars(distribution_vector(j)),'_s_',num2str(sample_vec(k))];
        figure('Name',fig_name)
        hold on;
        for t = 1:trials
            n = length(nse_lm{k,t});
            block_vec = linspace(1,n,n);
            plot(block_vec, nse_lm{k,t}, '-o')
        end
        xlabel('Block Number')
        ylabel('# of LM Per Block')
        bp = gca;
        if save_figs
            saveas(bp, fullfile(fig_dir, [fig_name, '.png']))
        end
    end

    % PDFs ----------------------------------------------------------------
    for k = 1:length(sample_vec)
        fig_name = [sub_dir,'_pdf_d_',convertStringsToChars(distribution_vector(j)),'_s_',num2str(sample_vec(k))];
        figure('Name',fig_name)
        hold on;
        for i = 1:trials
            nmem_h = plot(nmem_pdfs{j,i,k,1}, nmem_pdfs{j,i,k,2}, '-','Color',[0 0 0]+0.8, DisplayName='pdf nmem');
            nse_h = plot(nse_pdfs{j,i,k,1}, nse_pdfs{j,i,k,2}, '-k', DisplayName='pdf nap');
        end
        xlabel('x')
        ylabel('f(x)')
        if max(nse_pdfs{j,1,1,1}) > 20
            xlim([0,20])
        end
        if max(nse_pdfs{j,1,1,2}) > 1
            ylim([0,6])
        else
            ylim([0,1])
        end
        legend([nmem_h(1),nse_h(1)],'Location','best')
        bp = gca;
        if save_figs
            saveas(bp, fullfile(fig_dir, [fig_name, '.png']))
        end
    end
    
%     % CDFs ----------------------------------------------------------------
%     for k = 1:length(sample_vec)
%         figure('Name',['CDFs for sample size ',num2str(sample_vec(k))])
%         hold on;
%         for i = 1:trials
%             plot(nmem_cdfs{j,i,k,1}, nmem_cdfs{j,i,k,2}, '-b', DisplayName='pdf nmem')
%             plot(nse_cdfs{j,i,k,1}, nse_cdfs{j,i,k,2}, '-r', DisplayName='pdf nap')
%         end
%         xlabel('x')
%         ylabel('f(x)')
%     end
%  
%     % SQRs ----------------------------------------------------------------
%     for k = 1:length(sample_vec)
%         figure('Name',['SQRs for sample size ',num2str(sample_vec(k))])
%     
%         hold on;
%         smallN = 256;
%         smallN2 = 258;
%         graymax = 240;
%         range = 0:1/(smallN+1):1;
%         muLD = range*(smallN + 1) / (smallN + 1);
%         lemonDrop = sqrt(muLD.*(1-muLD)) * 3.4;
%         sampleCount2 = (smallN + 2):-1:1;
%         colorRange = (255-graymax)*sampleCount2/(smallN + 2);
%         base = repmat(graymax, smallN + 2, 1);
%         col = (base + colorRange') / 255;
%         rgb = [col col col];
%         count2 = 1;
%         for ii = ceil(smallN2/2):smallN2-1
%             ix = [ii ii+1 smallN2-ii smallN2-ii+1];
%             fill(range(ix), lemonDrop(ix), rgb(count2, :),'edgecolor','none')
%             fill(range(ix), -lemonDrop(ix), rgb(count2, :),'edgecolor','none')
%             count2 = count2 + 2;
%         end
%         plot(muLD,lemonDrop,'k--');
%         plot(muLD,-lemonDrop,'k--');
%         %             h = plot(nap(d_indx,n_indx).u, nap(d_indx,n_indx).sqr, 'DisplayName', '$\hat{sqr}_i(x)$');
%    
%         for i = 1:trials
%             nse_h = plot(nmem_sqrs{j,i,k,1}, nmem_sqrs{j,i,k,2},'-r', 'DisplayName', '$\hat{sqr}_i^{NAP}(x)$');
%             nmem_h = plot(nse_sqrs{j,i,k,1}, nse_sqrs{j,i,k,2},'-b', 'DisplayName', '$\hat{sqr}_i^{NMEM}(x)$');
%         end
%         
%         bp = gca;
%         xlabel('$x$','Interpreter','latex')
%         ylabel('$\hat{SQR}(x)$','Interpreter','latex')
%         %             legend([nse_h(1),nmem_h(1), g],'Interpreter','latex')
%         xlabel('$x$','Interpreter','latex')
%         ylabel('$sqr(x)$','Interpreter','latex')
%     end





end






