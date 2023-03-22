clc; clear; close all;

addpath("functions/")
addpath("functions_plotting/")



%%%%%%%%%%%%%%%%%%%%%%%%%%%%% To Plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

plot_mse_dist = false;
table_name = 'mse_kl_1.dat';
write_dir = fullfile('data_large_n','kl_mse_cpu_40_t_1');
status = mkdir(write_dir);

plot_pdf = true;
plot_cdf = false;
plot_sqr = true;
plot_heavy = false;
save_figs = false;

calc_q = false;

% figure directory
fig_dir = fullfile('figures_manuscript_large_n','1_trials');
status = mkdir(fig_dir);

% choose to visualise figures or not
fig_plot = 'on';

%%%%%%%%%%%%%%%%%%%%%%%%%%% Import data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Create actual distrobution onjects
actual = distributions;

% Path to directory
% dir_name = fullfile('data','pdf_estimates');
dir_name = fullfile('data_large_n','cpu_40_t_1');
write_name = fullfile('data_large_n','kl_mse_cpu_40_t_1');
status = mkdir(write_name);

% Define the etimates to plot

% Sample range
% n_vec = 2.^[8,9,10,11,12,13,14,15,16,17,18];
n_vec = 2.^[14];
n_vec = 2.^[8,9,10,11,12,13,14,15,16,17,18,19,20,21,22];
n_vec = 2.^[8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27];
% n_vec = 2.^[10,11,12,13,14,15,16,17,18,19,20,21,22];
% n_vec = 2.^[8,,11,12,13];
% n_vec = 2.^[10,15];
% Distribution range
% d_vec = ["Beta-a0p5-b1p5","Beta-a2-b0p5","Beta-a0p5-b0p5","Generalized-Pareto","Stable"];
% d_vec = ["Trimodal-Normal","Uniform","Normal","Uniform-Mix","Beta-a0p5-b1p5","Beta-a2-b0p5","Beta-a0p5-b0p5","Generalized-Pareto","Stable"];
d_vec = ["Trimodal-Normal","Uniform","Normal","Beta-a0p5-b1p5","Beta-a2-b0p5","Beta-a0p5-b0p5","Generalized-Pareto","Stable"];
d_vec = ["Beta-a0p5-b1p5"];
% d_vec = ["Trimodal-Normal","Uniform","Normal"];
% d_vec = ["Beta-a0p5-b1p5","Beta-a2-b0p5","Beta-a0p5-b0p5","Generalized-Pareto"];
% d_vec = ["Normal"];

% Define labels for figures
% names = {'Uniform-Mix', 'Generalized-Pareto', 'Stable',...
%     'Tri-Modal-Normal', 'Normal', 'Uniform', 'Beta(0.5,1.5)',...
%     'Beta(2,0.5)', 'Beta(0.5,0.5)'}';
% names = {'Beta(0.5,1.5)','Beta(2,0.5)', 'Beta(0.5,0.5)', 'Generalized-Pareto', 'Stable'}';
names = ["Trimodal-Normal","Uniform","Normal","Beta(0.5,1.5)","Beta(2,0.5)", "Beta(0.5,0.5)", "Generalized-Pareto","Stable"];
names = ["Beta(0.5,1.5)"];

% d_vec = ["Beta-a0p5-b0p5","Generalized-Pareto"];
% names = ["Beta(0.5,0.5)", "Generalized-Pareto"];
%
% d_vec = ["Uniform","Trimodal-Normal","Beta-a0p5-b1p5","Beta-a2-b0p5","Beta-a0p5-b0p5"];
% names = ["Uniform","Trimodal-Normal","Beta(0.5,1.5)","Beta(2,0.5)", "Beta(0.5,0.5)"];
%
% d_vec = ["Beta-a0p5-b1p5","Beta-a2-b0p5","Beta-a0p5-b0p5"];
% names = ["Beta(0.5,1.5)","Beta(2,0.5)", "Beta(0.5,0.5)"];

% Trials per sample
trials = 1;
plt_trials = 1;
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

            % x
            nse_data = nse_pdf_data;

            nap_x = nse_data{2,k}'; % make sure column vector

            % pdf
            nap_pdf = nse_data{3,k}'; % make sure column vector

            % cdf
            nap_cdf = utils_analysis.trapz(nse_data{2,k}, nse_data{3,k});

            % calcualte kl,mse,quantile -----------------------------------

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

            filename = sprintf('mse_per_q_%s_n_%s_t_%s.mat',d_vec(i), num2str(n_vec(j)), num2str(trials));

            save(fullfile(write_dir,['nse_', filename]), 'mse_dist_nse_per_q')
        end


    end

    % build table


    % MSE
    temp = utils.reshape_groups(n_vec',mse_dist_nse);

    sample_power = temp(:,1);
    mse_total = temp(:,2);

    % KL
    temp = utils.reshape_groups(n_vec',kl_dist_nse);

    kl_total = temp(:,2);

    % Distributions
    distribution = repelem(d_vec(i), length(temp(:,2)))';
    name = repelem(names(i), length(temp(:,2)))';

    nse_label = repelem(["NAP"], size(utils.reshape_groups(n_vec',mse_dist_nse), 1));


    estimator = nse_label';


    dist_table = table(distribution, name, estimator, sample_power, mse_total,kl_total);
    dist_table = splitvars(dist_table);

    % append table per distribution to global table containing data
    % for all distributions
    global_table = [global_table; dist_table];
end


writetable(global_table,fullfile(write_dir,table_name))

test = 'test'
