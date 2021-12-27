% use when driver.m is a script and not a function
clc;clear;
% class assignment
actual = distributions;
actual.generate_data = false;

% SUBSAMPLING PARAMETERS---------------------------------------------------

BootSample = "off";
if BootSample == "off"
    % percentage of sample used to create subsample
    percSample = 1;
    % number of subsamples to generate
    numSubs = 1;
else
    % percentage of sample used to create subsample
    percSample = 0.4;
    % number of subsamples to generate
    numSubs = 30;
end

%--------------------------------------------------------------------------

tic
% User Options ============================================================
% script switching board
estimator_call_flag =       true;   %<- true/false call SE on/off
estimator_plot_flag =       false;   %<- true/false plot SE results on/off
data_type_flag =            true;   %<- true/false integer powers of 2/real powers of 2
save_graphics =             false;   %<- true/false save .png of plots on/off
% rndom data generation parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
max_pow =                   10; %<---- maximum exponent to generate samples
min_pow =                   8; %<---- minimum exponent to generate samples
trials =                    4;  %<--- trials to run to generate heuristics for programs
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
distribution_vector = ["Generalized-Pareto","Stable","Trimodal-Normal","Normal", "Uniform","Beta-a0p5-b1p5","Beta-a2-b0p5","Beta-a0p5-b0p5"];
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
% NOTE: only works for division with no remainders
fail_nmem = zeros(trials,(max_pow-min_pow)/step,length(distribution_vector));
fail_se = zeros(trials,(max_pow-min_pow)/step,length(distribution_vector));

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
    T_track = zeros(max_pow-min_pow,1);
    % track BR per branch level
    BR_track = zeros(max_pow-min_pow,trials);
    % track BR at the zeroth  branch level
    BR0_track = zeros(max_pow-min_pow,trials);
    %
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
    dist_BR0_track = cell(length(distribution_vector));
    length(distribution_vector)
    
    disp(max_pow-min_pow)
    
    dist_BR_track = cell(max_pow-min_pow);
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

            % run estimator and store data
            if estimator_call_flag
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
                
                % while loop to replace non-finite values with finite values

                %                 pd = rndom.distInfo;
                %                 non_finite_vals = sum(~isfinite(sample));
                %                 % check sample
                %                 for inf_index = 1:non_finite_vals
                %                     index = isinf(sample);
                %                     num_inf = non_finite_vals;
                %
                %                     new_data = inf;
                %                     disp(['Old: ', num2str(sample(inf_index))])
                %                     while isinf(new_data)
                %                         new_data = random(pd);
                %                         sample(inf_index) = new_data;
                %                     end
                %                     disp(['New: ', num2str(sample(inf_index))])
                %                 end
                %

                % mixture model does not return random sample
                [fail_code,x,SE_pdf,SE_cdf,SE_u,SE_SQR,nBlocks,Blacklist,rndom.Ns,binrndom.Ns, max_LG, sum_LG,T,BRlevel,BR0]...
                    = stitch_pdf(sample,rndom.filename,send_file_name,actual.min_limit,actual.max_limit,p);

                fail_se(i,k,j) = fail_code;

                % track T,BR per trial

                T_track(k) = T;
                BR0_track(k,i) = BR0;
                BR_track(k,i) = sum(BRlevel{end,1});

                    % store subsample estimate data
                estimate_data.x = x;
                estimate_data.pdf = SE_pdf;
                estimate_data.u = SE_u;
                estimate_data.sqr = SE_SQR;


                tcpuSE = cputime-tintialSE;
                tintialNMEM = cputime;

                %-- NMEM start
                %[failed, x_NMEM, pdf_NMEM, cdf_NMEM,sqr_NMEM, ~,lagrange] = EstimatePDF(sample);
                try
                    [failed, x_NMEM, pdf_NMEM, cdf_NMEM,sqr_NMEM, ~,lagrange] = EstimatePDF(sample);
                    fail_nmem(i,k,j) = 0;
                    dlmwrite(['NMEM_pdf_',rndom.filename,'.dat'],[x_NMEM, pdf_NMEM],'Precision',12)
                    dlmwrite(['NMEM_cdf_',rndom.filename,'.dat'],[x_NMEM, cdf_NMEM],'Precision',12)

                    n = length(sqr_NMEM);
                    dx = 1 / (n + 1);
                    u_NMEM = dx:dx:(n * dx);

                    dlmwrite(['NMEM_sqr_',rndom.filename,'.dat'],[u_NMEM', sqr_NMEM],'Precision',12)
                catch
                    warning('Problem using function.  Assigning a value of 0.');
                    lagrange = 0;
                    x_NMEM = linspace(min(sample),max(sample),length(sample))';
                    pdf_NMEM = 0*ones(length(sample),1);
                    cdf_NMEM = 0*ones(length(sample),1);
                    fail_nmem(i,k,j) = 1;
                end
                %----------------------------------------------------------
                tcpuNMEM = cputime-tintialNMEM ;
                disp(['NMEM elapsed time is ',num2str(tcpuNMEM),' seconds.'])
                %==========================================================
                % % % % % % % % end of estimate % % % % % % % % %
                %==========================================================

                % calculate LG multipliers
                max_LG_se(k,i) = sum(max_LG);
                max_LG_nmem(k,i) = sum(length(lagrange));

                % calculate kl values
                kl_info_se = actual;
                kl_info_se.x = x;
                kl_info_se = dist_list(kl_info_se);

                kl_info_nmem = actual;
                kl_info_nmem.x = x_NMEM';
                kl_info_nmem = dist_list(kl_info_nmem);
                kl_info_nmem;
                %----------------------------------------------------------

                SE_test_x = interp1(kl_info_se.x,kl_info_se.pdf_y,x);
                NMEM_test_x = interp1(kl_info_nmem.x,kl_info_nmem.pdf_y,x_NMEM');

                if sum(~isfinite(interp1(kl_info_se.x,kl_info_se.pdf_y,x)))...
                        || sum(~isfinite(SE_pdf))... % for some reason indefinite values
                        || sum(~isfinite(pdf_NMEM))...
                        || sum(~isfinite(interp1(kl_info_nmem.x,kl_info_nmem.pdf_y,x_NMEM')))
                    disp('non-finite values')
                end
                actual_se_pdf = interp1(kl_info_se.x,kl_info_se.pdf_y,x);
                actual_nmem_pdf = interp1(kl_info_nmem.x,kl_info_nmem.pdf_y,x_NMEM');

                % find rand sample to find kl for
                [s_se,idx1] = datasample(SE_pdf,length(x_NMEM));
                [s_nmem,idx2] = datasample(pdf_NMEM',k);

                if max(x_NMEM) > max(x) && min(x) > min(x_NMEM)
                    kl_se_test = sum(SE_pdf' - interp1(actual.x',actual.pdf_y',x'));
                    kl_nmem_test = sum(interp1(x_NMEM,pdf_NMEM,x') - interp1(actual.x',actual.pdf_y',x'));
                elseif max(x_NMEM) < max(x) && min(x) < min(x_NMEM)
                    kl_se_test = sum(interp1(x',SE_pdf',x_NMEM) - interp1(actual.x',actual.pdf_y',x_NMEM));
                    kl_nmem_test = sum(pdf_NMEM - interp1(actual.x',actual.pdf_y',x_NMEM));
                end

                % read in actual distribution

                rndom.dist_list();

                % Create final answer file
                rndom.filename = sprintf(['D_', char(rndom.dist_name),'_T_','%d', '_S_','%d'],i, rndom.Ns);

                Sanswer = [x',SE_pdf'];
                S_sqr = [SE_u',SE_SQR'];
                S_cdf = [x',SE_cdf'];
                if fail_code == 0
                    dlmwrite(['SE_pdf_',rndom.filename,'.dat'],Sanswer, 'delimiter',' ','precision',12)
                    dlmwrite(['SE_sqr_',rndom.filename,'.dat'],S_sqr, 'delimiter',' ','precision',12)
                    dlmwrite(['SE_cdf_',rndom.filename,'.dat'],S_cdf, 'delimiter',' ','precision',12)
                end
            end

            % SE PLOTS-----------------------------------------------------
            % set figure styles
            publicationQuality();

            % SQR plot function
            [u_estimate,sqr_estimate] = misc_functions.sqr(x,SE_pdf,sample);
            [u_NMEM,sqr_NMEM] = misc_functions.sqr(x_NMEM,pdf_NMEM,sample);

            %{
            if estimator_plot_flag
                figure('Name','Bootstrap Estimates')
                hold on
                plot(x,SE_pdf,'Color',[0,0,0],'DisplayName','Average');
                plot(x_NMEM,pdf_NMEM,'Color',[1,0,0],'DisplayName','Average');
                plot(actual.x,actual.pdf_y,'--k')
                if max(SE_pdf) > 2
                    ylim([0,6])
                else
                    ylim([0,1])
                end
                xlim([actual.min_limit,actual.max_limit])
                ylabel('$PDF$','Interpreter','latex')
                xlabel('x','Interpreter','latex')
                if save_graphics
                    png_file = strcat('boot_PDF_D_',char(actual.dist_name),'T_',int2str(i),'S_',int2str(sample_vec(k)),'_',int2str(percSample),'_',int2str(numSubs),'.png');
                    saveas(gcf,png_file)
                    fig_file = strcat(png_file,'.eps');
                    print(fig_file,'-depsc')
                end

                %sqr ------------------------------------------------------
                figure('Name',['SE SQR rndom.Ns: ',...
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
                ldg2(1) = plot(u_estimate,sqr_estimate,'-k','DisplayName','BSSE');
                ldg2(2) = plot(u_NMEM,sqr_NMEM,'-r','DisplayName','BSSE');
                xlim([0,1])
                ylabel('$SQR$','Interpreter','latex')
                xlabel('u','Interpreter','latex')
                if save_graphics
                    graphics_name = ['SQR_D_',...
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



            end
            %}

            % heavy tail analysis -------------------------------------
            %{
            basex = 1;%log(exp(1));
            basey = log(1000);
            upperExtend = 1;
            lowerExtend = 1;

            posXnmem = x_NMEM(x_NMEM>=0);
            posXse = x(x>=0);
            posCDFnmem = cdf_NMEM(x_NMEM>=0);
            posCDFse = SE_cdf(x>=0);

            xHTnmem = log(posXnmem)/basex;
            xHTse = log(posXse)/basex;
            yHTnmem = log(1-posCDFnmem)/basey;
            yHTse = log(1-posCDFse)/basey;

            HTactual = actual;

            if max(posXnmem) > max(posXse) && min(posXse) > min(posXnmem)
                HTactual.x = posXnmem;
                HTactual = dist_list(HTactual);
            elseif max(posXnmem) < max(posXse) && min(posXse) < min(posXnmem)
                HTactual.x = posXse;
                HTactual = dist_list(HTactual);
            elseif max(posXnmem) < max(posXse) && min(posXse) > min(posXnmem)
                HTactual.x = linspace(min(posXnmem),max(posXse),10000);
                HTactual = dist_list(HTactual);
            elseif max(posXnmem) > max(posXse) && min(posXse) < min(posXnmem)
                HTactual.x = linspace(min(posXse),max(posXnmem),10000);
                HTactual = dist_list(HTactual);
            end

            xHTact = log(HTactual.x)/basex;

            yHTact = log(1-HTactual.cdf_y)/basey;

            figure('Name',['heavy tail compare',char(HTactual.dist_name),...
                'T_',int2str(trials),...
                'S_',int2str(sample_vec(k)),...
                '_',int2str(percSample),...
                '_',int2str(numSubs)])
            hold on
            plot(xHTnmem,yHTnmem,'.r')
            plot(xHTse,yHTse,'.b')
            plot(xHTact,yHTact,'-k')

            xlim([lowerExtend*min(xHTact),upperExtend*max(xHTact)]);
            ylim([lowerExtend*min(yHTact),upperExtend*max(yHTact)]);
            ylabel('$\frac{ln(1 - \hat{F}(x))}{3ln(10)}$','Interpreter','latex')
            xlabel('$ln(x)$','Interpreter','latex')
            legend('$\hat{f}_{NMEM}(x)$','$\hat{f}_{SE}(x)$','Interpreter','latex','Location','SouthWest')
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

            % store pdf/sqr solutions for all dist/trials/samples
            all_se{i,k}.x(:,1) = x;
            all_se{i,k}.pdf(:,1) = SE_pdf;
            all_nmem{i,k}.x(:,1) = x_NMEM;
            all_nmem{i,k}.pdf(:,1) = pdf_NMEM;
            all_se{i,k}.u(:,1) = u_estimate;
            all_se{i,k}.sqr(:,1) = sqr_estimate;
            all_nmem{i,k}.u(:,1) = u_NMEM;
            all_nmem{i,k}.sqr(:,1) = sqr_NMEM;

            disp([char(actual.dist_name),...
                ', Trial: ',num2str(i),'/', num2str(trials), ...
                ' sample size: ',num2str(sample_vec(k)), ...
                ' failSE: ', num2str(fail_se(i,k,j)), ' failNMEM: ', num2str(fail_nmem(i,k,j))])

            % store time of computation
            cpu_vec_se(k,i) = tcpuSE;
            cpu_vec_nmem(k,i) = tcpuNMEM;

            %{
            % SQR bootstrap plot function
            figure('Name',['Boot_brndom.Ns: ',num2str(bootrandomNs),' rndom.Ns: ', num2str(rndom.Ns),' blocks: ',num2str(nBlocks),' numSubs: ',num2str(numSubs),' percSample: ',num2str(percSample)])
            hold on
            l = zeros(1,2);
            %l(1) = plot(estimate_data(1).x,estimate_data(1).pdf,'Color',[.8,.8,.8],'DisplayName','sub-samples');
            plot(estimate_data(1).x,estimate_data(1).pdf,'Color',[.8,.8,.8],'DisplayName','sub-samples');
            for s = 2:size(estimate_data,2)
                plot(estimate_data(s).x,estimate_data(s).pdf,'Color',[.8,.8,.8])
            end
            l(2) = plot(x_NMEM,pdf_NMEM,'-r','DisplayName','NMEM');
            l(1) = plot(xinterp,AvgEst,'-b','DisplayName','SE');
            plot(actual.x,actual.pdf_y,'--k')
            if max(SE_pdf) > 2
                ylim([0,6])
            else
                ylim([0,1])
            end
            xlim([actual.min_limit,actual.max_limit])
            ylabel('$PDF$','Interpreter','latex')
            xlabel('x','Interpreter','latex')
            %legend(l(1))
            legend(l(1:2))
            if save_graphics
                png_file = strcat('boot_PDF_D_',char(actual.dist_name),'T_',int2str(i),'S_',int2str(sample_vec(k)),'_',int2str(percSample),'_',int2str(numSubs),'.png');
                saveas(gcf,png_file)
                fig_file = strcat(png_file,'.eps');
                print(fig_file,'-depsc')
            end
            %}
        end

        se_time_per_trial = horzcat(sample_vec',cpu_vec_se(:,i));
        nmem_time_per_trial = horzcat(sample_vec',cpu_vec_nmem(:,i));
        
        [p,q] = size(se_time_per_trial); 
        BoxCPUtimeSE(end-p+1:end, end-q+1:end) = nmem_time_per_trial;
        [p,q] = size(nmem_time_per_trial); 
        BoxCPUtimeNMEM(end-p+1:end, end-q+1:end) = nmem_time_per_trial;
        
        rndom.filename = sprintf(['D_', char(actual.dist_name),'_T_','%d'],i);
    end

    disp(j)
    disp(size(dist_BR0_track))
    disp(size(dist_BR0_track{j}))
    disp(size(BR0_track))
    
    dist_BR0_track{j} = BR0_track;
    dist_BR_track{j} = BR_track;

    %{
    figure('Name',['cpu Box Plot SE: ',char(actual.dist_name)])
    boxplot(BoxCPUtimeSE(:,2),BoxCPUtimeSE(:,1))
    xlabel('CPU Time','interpreter','latex')
    ylabel('Sample Size')
    if save_graphics
        graphics_name = ['BoxCPUse_',...
            char(actual.dist_name),...
            'T_',int2str(trials)];
        png_file = strcat(graphics_name,'.png');
        saveas(gcf,png_file)
        fig_file = strcat(graphics_name,'.eps');
        print(fig_file,'-depsc')
    end

    labels = {'$2^{15}$','$2^{17}$'};
%
%     figure('Name','SE MSE for estimates')
%     boxplot(MSE_se,sample_data_box,'Labels',labels)
%     bp = gca;
%     bp.XAxis.TickLabelInterpreter = 'latex';

    figure('Name',['cpu Box Plot NMEM: ',char(actual.dist_name)])
%     boxplot(BoxCPUtimeSE(:,2),BoxCPUtimeSE(:,1),'Labels',labels)
    boxplot(BoxCPUtimeSE(:,2),BoxCPUtimeSE(:,1))
    bp = gca;
    bp.XAxis.TickLabelInterpreter = 'latex';
    xlabel('CPU Time','interpreter','latex')
    ylabel('Sample Size')
    if save_graphics
        graphics_name = ['BoxCPUnmem_',...
            char(actual.dist_name),...
            'T_',int2str(trials)];
        png_file = strcat(graphics_name,'.png');
        saveas(gcf,png_file)
        fig_file = strcat(graphics_name,'.eps');
        print(fig_file,'-depsc')
    end
    %}

   %{
    for k = 1:length(sample_vec)
        % pdf
        figure('Name',['plot all estimates SE: ',char(actual.dist_name),'-',num2str(sample_vec(k))])
        hold on
        for i = 1:trials
            plot(all_se{i,k}.x(:,1), all_se{i,k}.pdf(:,1))
        end
        plot(actual.x,actual.pdf_y,'--k')
        if max(all_se{i,k}.pdf(:,1)) > 2
            ylim([0,6])
        else
            ylim([0,1])
        end
        xlim([actual.min_limit,actual.max_limit])
        if save_graphics
            graphics_name = ['pdf_all_SE_',...
                char(actual.dist_name),...
                'S_',num2str(sample_vec(k))];
            png_file = strcat(graphics_name,'.png');
            saveas(gcf,png_file)
            fig_file = strcat(graphics_name,'.eps');
            print(fig_file,'-depsc')
        end

        figure('Name',['plot all estimates NMEM: ',char(actual.dist_name),'-',num2str(sample_vec(k))])
        hold on
        for i = 1:trials
            plot(all_nmem{i,k}.x(:,1), all_nmem{i,k}.pdf(:,1))
        end
        plot(actual.x,actual.pdf_y,'--k')
        if max(all_nmem{i,k}.pdf(:,1)) > 2
            ylim([0,6])
        else
            ylim([0,1])
        end
        xlim([actual.min_limit,actual.max_limit])
        if save_graphics
            graphics_name = ['pdf_all_NMEM_',...
                char(actual.dist_name),...
                'S_',num2str(sample_vec(k))];
            png_file = strcat(graphics_name,'.png');
            saveas(gcf,png_file)
            fig_file = strcat(graphics_name,'.eps');
            print(fig_file,'-depsc')
        end

        % sqr
        figure('Name',['plot all SQR SE: ',char(actual.dist_name),'-',num2str(sample_vec(k))])
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
        for i = 1:trials
            plot(all_se{i,k}.u(:,1),all_se{i,k}.sqr(:,1),'DisplayName','SE');
        end
        xlim([0,1])
        ylabel('$SQR$','Interpreter','latex')
        xlabel('u','Interpreter','latex')
        if save_graphics
            graphics_name = ['plot_all_SQR_SE_D_',...
                char(actual.dist_name),...
                'T_',int2str(trials),...
                'S_',int2str(sample_vec(k)),...
                'perc_',int2str(percSample),...
                'sub_',int2str(numSubs)];
            png_file = strcat(graphics_name,'.png');
            saveas(gcf,png_file)
            fig_file = strcat(graphics_name,'.eps');
            print(fig_file,'-depsc')
        end

        figure('Name',['plot all SQR NMEM: ',char(actual.dist_name),'-',num2str(sample_vec(k))])
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
        for i = 1:trials
            plot(all_nmem{i,k}.u(:,1),all_nmem{i,k}.sqr(:,1),'DisplayName','SE');
        end
        xlim([0,1])
        ylabel('$SQR$','Interpreter','latex')
        xlabel('u','Interpreter','latex')
        if save_graphics
            graphics_name = ['plot_all_SQR_NMEM_D_',...
                char(actual.dist_name),...
                'T_',int2str(trials),...
                'S_',int2str(sample_vec(k)),...
                'perc_',int2str(percSample),...
                'sub_',int2str(numSubs)];
            png_file = strcat(graphics_name,'.png');
            saveas(gcf,png_file)
            fig_file = strcat(graphics_name,'.eps');
            print(fig_file,'-depsc')
        end

    end
    %}

    %--------------------------------------------------

    stdCPUtimeSE = horzcat(sample_vec',std(cpu_vec_se,[],2));
    stdCPUtimeNMEM = horzcat(sample_vec',std(cpu_vec_nmem,[],2));
    avgCPUtimeSE = horzcat(sample_vec',mean(cpu_vec_se,2));
    avgCPUtimeNMEM = horzcat(sample_vec',mean(cpu_vec_nmem,2));

    cpu.stdTimeSE{j} = horzcat(sample_vec',std(cpu_vec_se,[],2));
    cpu.stdTimeNMEM{j} = horzcat(sample_vec',std(cpu_vec_nmem,[],2));
    cpu.timeSE{j} = horzcat(sample_vec',mean(cpu_vec_se,2));
    cpu.timeNMEM{j} = horzcat(sample_vec',mean(cpu_vec_nmem,2));

    dlmwrite(['cpu_SE_',rndom.filename,'.dat'],horzcat(sample_vec',cpu_vec_se), 'delimiter',' ','precision',12)
    dlmwrite(['cpu_NMEM_',rndom.filename,'.dat'],horzcat(sample_vec',cpu_vec_nmem), 'delimiter',' ','precision',12)
    dlmwrite(['BRTrack_SE_',rndom.filename,'.dat'],horzcat(sample_vec',mean(dist_BR0_track{j},2)), 'delimiter',' ','precision',12)
    dlmwrite(['Failures_SE_',rndom.filename,'.dat'],horzcat(sample_vec',sum(fail_se(:,:,j),1)'/trials), 'delimiter',' ','precision',12)
    dlmwrite(['Failures_NMEM_',rndom.filename,'.dat'],horzcat(sample_vec',sum(fail_nmem(:,:,j),1)'/trials), 'delimiter',' ','precision',12)

end
toc
%{
figure('Name','testBR')
hold on
cc=lines(length(distribution_vector));
for o = 1:length(distribution_vector)
    h(o) = plot(sample_vec,mean(dist_BR0_track{o},2),'DisplayName',char(distribution_vector(o)),'Color',cc(o,:));
end
t(1) = plot(sample_vec,T_track,'--k','DisplayName','Threshold');
xlabel('2^x')
ylabel('Threshold')
legend([h,t])
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
%profile viewer
% END OF PROGRAM ==========================================================
