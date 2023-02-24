clc; clear; close all;

addpath("functions/")

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% To Plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

plot_mse_dist = true;
table_name = 'mse_kl_set0.dat';

plot_pdf = true;
plot_cdf = false;
plot_sqr = true;
plot_heavy = true;
save_figs = true;

% file name label
filename = 'mse_kl_set0.mat';

%%%%%%%%%%%%%%%%%%%%%%%%%%% Import data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Create actual distrobution onjects
actual = distributions;

% Path to directory
% dir_name = fullfile('data','pdf_estimates');
dir_name = fullfile('data','cpu_15_t_50_maxN_22_pdf_estimates');
save_dir = fullfile('data');

% Define the etimates to plot

% Sample range
% n_vec = 2.^[8,9,10,11,12,13,14,15,16,17,18];
n_vec = 2.^[8,9,10,11,12,13,14,15,16,17,18,19,20,21,22];
n_vec = 2.^[15,16];

% Distribution range
% d_vec = ["Trimodal-Normal","Uniform","Normal","Uniform-Mix","Beta-a0p5-b1p5","Beta-a2-b0p5","Beta-a0p5-b0p5","Generalized-Pareto","Stable"];
d_vec = ["Trimodal-Normal","Uniform","Normal"];

% Define labels for figures
% names = {'Uniform-Mix', 'Generalized-Pareto', 'Stable',...
%     'Tri-Modal-Normal', 'Normal', 'Uniform', 'Beta(0.5,1.5)',...
%     'Beta(2,0.5)', 'Beta(0.5,0.5)'}';
% names = {'Beta(0.5,1.5)','Beta(2,0.5)', 'Beta(0.5,0.5)', 'Generalized-Pareto', 'Stable'}';
names = ["Trimodal-Normal","Uniform","Normal"];

% Trials per sample
trials = 50;

% Initialize empty structure to store data in
n = length(n_vec);
d = length(d_vec);
t = length(d_vec);
nse = cell(d,n,t,7);
nmem = cell(d,n,t,7);

for i = 1:length(d_vec)
    for j = 1:length(n_vec)
        for k = 1:trials
                        
            % Make sure that distribution names are saved as character vectors

            % track calculations
            disp(['dist: ',num2str(i),'/',num2str(length(d_vec)),...
                ' s: ',num2str(j),'/',num2str(length(n_vec)),...
                ' t: ',num2str(k),'/',num2str(trials)])

            % distribution names
            nse{i,j,k,1} = d_vec(i);
            nmem{i,j,k,1} = d_vec(i);
            % sample size
            nse{i,j,k,2} = n_vec(j);
            nmem{i,j,k,2} = n_vec(j);

            % Load data from estimate directory
            load(fullfile(dir_name,['nse_',num2str(d_vec(i)),'_t_',num2str(trials),'_s_', num2str(n_vec(j)),'.mat']));
            load(fullfile(dir_name,['nmem_',num2str(d_vec(i)),'_t_',num2str(trials),'_s_', num2str(n_vec(j)),'.mat']));

            % x
            nse_data = nse_pdf_data;
            nmem_data = nmem_pdf_data;

            nse{i,j,k,3} = nse_data{2,k}'; % make sure column vector
            nmem{i,j,k,3} = nmem_data{2,k};

            % pdf
            nse{i,j,k,4} = nse_data{3,k}'; % make sure column vector
            nmem{i,j,k,4} = nmem_data{3,k};

            % cdf
            nse{i,j,k,5} = trapz(nse_data{2,k}, nse_data{3,k});
            nmem{i,j,k,5} = trapz(nmem_data{2,k}, nmem_data{3,k});

            % u and sqr
            [nse{i,j,k,6}, nse{i,j,k,7}] = utils.sqr(nse_data{2,k},nse_data{3,k},nse_data{1,k});
            [nmem{i,j,k,6}, nmem{i,j,k,7}] = utils.sqr(nmem_data{2,k},nmem_data{3,k},nmem_data{1,k});



        end
    end
end


save(fullfile(save_dir,['nse_', filename]), 'nse')
save(fullfile(save_dir,['nmem_', filename]), 'nmem')





