clc; clear; close all;

addpath("functions/")
addpath("functions_plotting/")

publicationQuality();


plot_mse_dist = false;
table_name = 'br0_table.dat';
data_dir = fullfile('data_cpu_20','theoretical_threshold');
fig_dir = fullfile('figures_manuscript','threshold');


status = mkdir(data_dir);
status = mkdir(fig_dir);

save_figs = false;
YexpScale = -2;

data = readtable(fullfile(data_dir,table_name));

col_names = data.Properties.VariableNames;

sample_vec = table2array(data(:,1));
distro_n = (length(col_names) - 1)/2;

names = ["Tri-Modal-Normal","Uniform", "Normal", "Beta(0.5,1.5)", "Beta(2,0.5)", "Beta(0.5,0.5)", "Generalized-Pareto", "stable"];
min_pow = min(sample_vec);
max_pow = max(sample_vec);

xT = linspace(2^min_pow,2^max_pow,10000);
% T = 2.6741.*xT.^(-0.626) + 0.005;
T = 2.*xT.^(0.33);

figure('Name','testBR')
hold on;
cc=lines(distro_n);
col_n = 2;
for j = 1:distro_n

    xorigin = sample_vec;
    yorigin = table2array(data(:,col_n));
    stdevorigin = table2array(data(:,col_n+1));
    column_name = data.Properties.VariableNames

    x = xorigin;
    y = yorigin;
    stdev = stdevorigin;

    curve1 = y + stdev;
    curve2 = y - stdev;

%     x2 = [x, fliplr(x)];
%     inBetween = [curve1, fliplr(curve2)];
%     f = fill(x2, inBetween,cc(j,:));
%     set(f,'facealpha',0.25)
%     set(f,'edgealpha',0)
    label = erase(char(col_names(col_n)),"Mean");

    h(j) = plot(x,y,'DisplayName',label,'Color',cc(j,:));

    col_n = col_n+2;
end
s = plot(log(xT)/log(2), T, 'k--','DisplayName','\Gamma');
bp = gca;
bp.XAxis.TickLabelInterpreter = 'latex';
% bp.xticklabel = labels;
bp.YAxis.Scale = "log";
xlabel('$log_{2}(x)$','Interpreter','latex')
ylabel('$\xi_0$','Interpreter','latex')
legend([s,h])



test = 'test'