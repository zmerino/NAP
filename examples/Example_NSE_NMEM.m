
% Clear unnecessary infromation
clc;clear all; close all;

% Add needed subdirectories
addpath("functions/")
addpath("compile_nmem/")

% class assignment used for general distribution plotting
actual = distributions;
actual.generate_data = false;

% User Options ============================================================

% script switching board %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data_type_flag =            true;   %<- true/false integer powers of 2/real powers of 2

% rndom data generation parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
max_pow =                   14; %<---- maximum exponent to generate samples
min_pow =                   14; %<---- minimum exponent to generate samples
trials =                    1   ;  %<--- trials to run to generate heuristics for programs
step =                      1;  %<---- control synthetic rndom samples to skip being created
temp_min_limit =            0; %<---- set upper limit for both
actual.min_limit =          temp_min_limit;  %<--- lower limit to plot
temp_max_limit =            1000; %<---- set upper limit for both
actual.max_limit =          temp_max_limit; %<--- upper limit to plot
x_resolution =              1000;
cpu_type =                   '\';%<--- '\' or '/' for windows or linux

% Example distribution to test %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
distribution_vector = ["Trimodal-Normal", "Normal"];

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

    % Create vector of  samples
    sample_vec = utils.sample_pow(min_pow,max_pow,data_type_flag,step);
    actual = actual.dist_list();
    rndom.Ns = sample_vec;

    % p-vector definition for Rtree function
    p = [1,0.55,1,0.33,2,ceil(0.0625*rndom.Ns^0.5),40];

    % Create rndom.filename for distribtuion
    rndom.filename = sprintf(['D_', char(actual.dist_name),'_T_','%d', '_S_','%d'],trials, rndom.Ns);

    % file path name
    send_file_name = ['D_',...
        char(actual.dist_name),...
        cpu_type,char(rndom.filename),...
        '.dat'];

    %==========================================================
    % % % % % % % start of estimates % % % % % % % % %
    %==========================================================

    %---------------------------------------------------------- NSE start
    sendFileName1 = ['D_',char(actual.dist_name),cpu_type,char(rndom.filename),'.dat'];
    rndom.randomVSactual = "random";
    rndom = dist_list(rndom);
    sample = rndom.rndData;

    % Use the following code to visualise internal workings of NSE()

    %     show_block_layout=true;
    %     visualize_interpolation=true;
    %     blockPdfs=true;
    %     savePNG=true;

    %     vis_list = num2cell([show_block_layout,visualize_interpolation,blockPdfs,savePNG]);

    %     stitchPdf = NSE(sample,rndom.filename, actual.min_limit,actual.max_limit,p,vis_list{:});

    % nse object instantiation
    nse = NSE;
    % Otherwise use defualt values and use the NSE class

    %     stitchPdf = nse.stitch(sample);
    stitchPdf = nse.stitch(sample);

    %---------------------------------------------------------- NMEM start
    try
        tintialNMEM = cputime;
        [failed, x_NMEM, pdf_NMEM, cdf_NMEM,sqr_NMEM, ~,lagrange_multipler] = EstimatePDF(sample);

        n = length(sqr_NMEM);
        dx = 1 / (n + 1);
        u_NMEM = dx:dx:(n * dx);
    catch
        warning('Problem using function.  Assigning a value of 0.');
        lagrange_multipler = 0;
        x_NMEM = linspace(min(sample),max(sample),length(sample))';
        pdf_NMEM = 0*ones(length(sample),1);
        cdf_NMEM = 0*ones(length(sample),1);
    end

    %==========================================================
    % % % % % % % % end of estimate % % % % % % % % %
    %==========================================================

    figure('Name','NSE and NMEM PDF')
    hold on;
    a = plot(actual.x, actual.pdf_y, 'k', 'DisplayName', 'Actual');
    b = plot(stitchPdf.sx, stitchPdf.sPDF, 'r', 'DisplayName', 'NSE');
    c = plot(x_NMEM, pdf_NMEM, 'b' , 'DisplayName', 'NMEM');
    bp = gca;
    bp.XAxis.TickLabelInterpreter = 'latex';
    xlabel('$log_{2}(N)$','Interpreter','latex')
    ylabel('$f(N)$','Interpreter','latex')
    legend()

end

