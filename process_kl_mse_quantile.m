clc; clear; close all;

addpath("functions/")
addpath("functions_plotting/")



%%%%%%%%%%%%%%%%%%%%%%%%%%%%% To Plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

plot_mse_dist = false;
table_name = 'mse_kl_100.dat';
write_dir = fullfile('data_4','kl_mse_data');
status = mkdir(write_dir);

plot_pdf = true;
plot_cdf = false;
plot_sqr = true;
plot_heavy = false;
save_figs = false;

calc_q = false;

% figure directory
fig_dir = fullfile('figures_manuscript','100_trials_v1');
status = mkdir(fig_dir);

% choose to visualise figures or not
fig_plot = 'on';

%%%%%%%%%%%%%%%%%%%%%%%%%%% Import data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Create actual distrobution onjects
actual = distributions;

% Path to directory
% dir_name = fullfile('data','pdf_estimates');
dir_name = fullfile('data_4','cpu_20_t_100');
write_name = fullfile('data_4','kl_mse_cpu_20_t_100');

% Define the etimates to plot

% Sample range
n_vec = 2.^[8,9,10,11,12,13,14,15,16,17,18,19,20,21,22];

d_vec = ["Trimodal-Normal","Uniform","Normal","Beta-a0p5-b1p5","Beta-a2-b0p5","Beta-a0p5-b0p5","Generalized-Pareto","Stable"];
names = ["Trimodal-Normal","Uniform","Normal","Beta(0.5,1.5)","Beta(2,0.5)", "Beta(0.5,0.5)", "Generalized-Pareto","Stable"];

% Trials per sample
trials = 100;
plt_trials = 10;
% number of quantile segments
qnum = 50;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% quantiles %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

actual = distributions;
actual.generate_data = false;

% Initialize empty structure to store data in
n = length(n_vec);
d = length(d_vec);
t = length(d_vec);
nap = cell(d,n,t,7);
nmem = cell(d,n,t,7);

mse_dist_nse_per_q = zeros(qnum,plt_trials);
mse_dist_nmem_per_q = zeros(qnum,plt_trials);
kl_dist_nmem_per_q = zeros(qnum,plt_trials);
kl_dist_nse_per_q = zeros(qnum,plt_trials);

global_table = table();

for i = 1:length(d_vec)
    % save distribution name to object
    actual.dist_name = d_vec(i);

    for j = 1:length(n_vec)
        for k = 1:trials

            % Make sure that distribution names are saved as character vectors

            % track calculations
            disp(['dist: ',num2str(i),'/',num2str(length(d_vec)),...
                ' s: ',num2str(j),'/',num2str(length(n_vec)),...
                ' t: ',num2str(k),'/',num2str(trials)])

            % Load data from estimate directory
            load(fullfile(dir_name,['nse_',num2str(d_vec(i)),'_t_',num2str(trials),'_s_', num2str(n_vec(j)),'.mat']));
            load(fullfile(dir_name,['nmem_',num2str(d_vec(i)),'_t_',num2str(trials),'_s_', num2str(n_vec(j)),'.mat']));

            % x
            nse_data = nse_pdf_data;
            nmem_data = nmem_pdf_data;

            nap_x = nse_data{2,k}'; % make sure column vector
            nmem_x = nmem_data{2,k};

            % pdf
            nap_pdf = nse_data{3,k}'; % make sure column vector
            nmem_pdf = nmem_data{3,k};

            % cdf
            nap_cdf = utils_analysis.trapz(nse_data{2,k}, nse_data{3,k});
            nmem_cdf = utils_analysis.trapz(nmem_data{2,k}, nmem_data{3,k});

            % calcualte kl,mse,quantile -----------------------------------
            % NMEM --------------------------------------------------------
            xs = nmem_x;
            fs = nmem_pdf;
            Fs = nmem_cdf;
            % get actual distribution
            actual.min_limit = min(xs);
            actual.max_limit = max(xs);
            actual.x = xs;
            actual = actual.dist_list();

            if sum(~isfinite(xs)) || sum(~isfinite(fs))|| sum(~isfinite(Fs))
                warning('non-finite values')
                test = 'test';
            end

            % MSE per quantile
            if calc_q
                [mse_dist_nmem_per_q(:,k), quantiles,xq,fq] = utils_analysis.quant_metric(xs, fs, Fs, qnum, d_vec(i),'MSE');
            end
            % MSE per trial
            mse_dist_nmem(j, k) = utils_analysis.mse(actual.pdf_y, fs);
            % KL per quantile
            if calc_q
                [kl_dist_nmem_per_q(:,k), ~,~,~] = utils_analysis.quant_metric(xs, fs, Fs, qnum, d_vec(i),'KL');
            end
            % KL per trial
            kl_dist_nmem(j, k) = utils_analysis.kl(actual.pdf_y, fs');


            % NAP ------------------------------------------------------
            %
            xs = nap_x;
            fs = nap_pdf;
            Fs = nap_cdf;
            % get actual distribution
            actual.min_limit = min(xs);
            actual.max_limit = max(xs);
            actual.x = xs;
            actual = actual.dist_list();

            % MSE per quantile
            if calc_q
                [mse_dist_nse_per_q(:,k),~,~,~] = utils_analysis.quant_metric(xs, fs, Fs, qnum, d_vec(i),'MSE');
            end
            % MSE per trial
            mse_dist_nse(j, k) = utils_analysis.mse(actual.pdf_y, fs);
            % KL per quantile
            if calc_q
                [kl_dist_nse_per_q(:,k), ~,~,~] = utils_analysis.quant_metric(xs, fs, Fs, qnum, d_vec(i),'KL');
            end
            % KL per trial
            % TODO: make sure sample data is transpose for single shot NAP
            try
                test = utils_analysis.kl(actual.pdf_y, fs);
                kl_dist_nse(j, k) = utils_analysis.kl(actual.pdf_y, fs);
            catch
                test = utils_analysis.kl(actual.pdf_y, fs');
                kl_dist_nse(j, k) = utils_analysis.kl(actual.pdf_y, fs');
            end

        end

        if calc_q
            filename = sprintf('kl_per_q_%s_n_%s_t_%s.mat',d_vec(i), num2str(n_vec(j)), num2str(trials));

            save(fullfile(write_dir,['nse_', filename]), 'kl_dist_nse_per_q')
            save(fullfile(write_dir,['nmem_', filename]), 'kl_dist_nmem_per_q')

            filename = sprintf('mse_per_q_%s_n_%s_t_%s.mat',d_vec(i), num2str(n_vec(j)), num2str(trials));

            save(fullfile(write_dir,['nse_', filename]), 'mse_dist_nse_per_q')
            save(fullfile(write_dir,['nmem_', filename]), 'mse_dist_nmem_per_q')
        end


    end

    % build table


    % MSE
    temp = vertcat(utils.reshape_groups(n_vec',mse_dist_nse),...
        utils.reshape_groups(n_vec',mse_dist_nmem));

    sample_power = temp(:,1);
    mse_total = temp(:,2);

    % KL
    temp = vertcat(utils.reshape_groups(n_vec',kl_dist_nse),...
        utils.reshape_groups(n_vec',kl_dist_nmem));

    kl_total = temp(:,2);

    % Distributions
    distribution = repelem(d_vec(i), length(temp(:,2)))';
    name = repelem(names(i), length(temp(:,2)))';

    nse_label = repelem(["NAP"], size(utils.reshape_groups(n_vec',mse_dist_nse), 1));
    nmem_label = repelem(["NMEM"], size(utils.reshape_groups(n_vec',mse_dist_nmem), 1));


    estimator = vertcat(nse_label', nmem_label');


    dist_table = table(distribution, name, estimator, sample_power, mse_total,kl_total);
    dist_table = splitvars(dist_table);

    % append table per distribution to global table containing data
    % for all distributions
    global_table = [global_table; dist_table];
end


writetable(global_table,fullfile(write_dir,table_name))
