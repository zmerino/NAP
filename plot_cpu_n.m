
close all;

addpath("functions/")
addpath("functions_plotting/")

publicationQuality();
save_figs = true;

% dir_name = fullfile('figures_manuscript','cpu_v1');
dir_name = fullfile('figures_manuscript_final_test','run_times');
status = mkdir(dir_name);

% Import data table
filename = fullfile('data_large_n', 'meta_data_1.dat');
data = readtable(filename);

% temporary fix: added inpuut sample size to data table
table_size = size(data ); 
rows = table_size(1);
idx = 1;
for i = 1:rows
    if mod(i+6, 27) == 0
        idx = 1;
    end
    data.sample_power(i) = 2^(idx+7);
    idx = idx + 1;
end

% Define labels for figures
names = ["Beta(0.5,1.5)", "Generalized-Pareto"];
labels = {'$2^{10}$', '$2^{11}$', '$2^{12}$',...
    '$2^{13}$', '$2^{14}$', '$2^{15}$', '$2^{16}$', '$2^{17}$',...
    '$2^{18}$', '$2^{19}$', '$2^{20}$', '$2^{21}$', '$2^{22}$',...
    '$2^{23}$', '$2^{24}$', '$2^{25}$', '$2^{26}$', '$2^{27}$'};
label_val = 2.^[8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27];


fig_name = 'CPU Time';
figure('Name',fig_name)
b = boxchart(log(data.sample_power)/log(2), data.cpu_time, 'GroupByColor',data.estimator);
bp = gca;
bp.XAxis.TickLabelInterpreter = 'latex';
xlabel('$log_{2}(N)$','Interpreter','latex')
ylabel('CPU Time','Interpreter','latex')
legend('Location','northwest')
if save_figs
    saveas(bp, fullfile('figures', [fig_name, '.png']))
end


fig_name = 'Wall Time';
figure('Name',fig_name)
b = boxchart(log(data.sample_power)/log(2), data.wall_time, 'GroupByColor',data.estimator);
bp = gca;
bp.XAxis.TickLabelInterpreter = 'latex';
xlabel('$log_{2}(N)$','Interpreter','latex')
ylabel('Wall Time','Interpreter','latex')
legend('Location','northwest')
if save_figs
    saveas(bp, fullfile('figures', [fig_name, '.png']))
end

fig_name = 'Failure Rate';
figure('Name',fig_name)
b = boxchart(log(data.sample_power)/log(2), data.fail, 'GroupByColor',data.estimator);
bp = gca;
bp.XAxis.TickLabelInterpreter = 'latex';
xlabel('$log_{2}(N)$','Interpreter','latex')
ylabel('Failure Rate','Interpreter','latex')
legend('Location','northwest')
if save_figs
    saveas(bp, fullfile(dir_name, [fig_name, '.png']))
end

% line plots ==============================================================

xlimits = [7,28];
ylimits = [-4,12];

fig_name = 'line_NAP Wall Time per Distribution';
figure('Name',fig_name)
% compute average ---------------------------------------------------------
data_wall = convertvars(data,["distribution","name","estimator","sample_power"],"categorical");
mean_data__parallel = groupsummary(data_wall,["distribution","name","estimator","sample_power"],"mean",["wall_time","cpu_time"]);
T = mean_data__parallel;
% Group based on first column - gender column
% Modify 1 to appropriate column number
G = findgroups(T{:, 1});
% Split table based on first column - gender column
T_split = splitapply( @(varargin) varargin, T , G);
% Allocate empty cell array fo sizxe equal to number of rows in T_Split
subTables = cell(size(T_split, 1));
% Create sub tables
for i = 1:size(T_split, 1)
subTables{i} = table(T_split{i, :}, 'VariableNames', ...
T.Properties.VariableNames);
end

hold on;
for i = 1:length(subTables)
    data_per_dist = subTables{i,1};
    x = double(string(data_per_dist.sample_power));
    y = double(data_per_dist.mean_wall_time);    

    x = log(x)/log(2);
    y = log(y);

    plot(x,y,'.-','DisplayName',string(data_per_dist.name(1)))
end
bp = gca;
bp.XAxis.TickLabelInterpreter = 'latex';
xlabel('$log_{2}(N)$','Interpreter','latex')
ylabel('log(CPU) Time','Interpreter','latex')
legend('Location','northwest')
xlim(xlimits)
ylim(ylimits)
if save_figs
    saveas(bp, fullfile(dir_name, [fig_name, '.png']))
end


xlimits = [7,28];
ylimits = [-4,12];

fig_name = 'line_NAP CPU per Distribution';
figure('Name',fig_name)
% compute average ---------------------------------------------------------
data_cpu = convertvars(data,["distribution","name","estimator","sample_power"],"categorical");
mean_data__parallel = groupsummary(data_cpu,["distribution","name","estimator","sample_power"],"mean",["wall_time","cpu_time"]);
T = mean_data__parallel;
% Group based on first column - gender column
% Modify 1 to appropriate column number
G = findgroups(T{:, 1});
% Split table based on first column - gender column
T_split = splitapply( @(varargin) varargin, T , G);
% Allocate empty cell array fo sizxe equal to number of rows in T_Split
subTables = cell(size(T_split, 1));
% Create sub tables
for i = 1:size(T_split, 1)
subTables{i} = table(T_split{i, :}, 'VariableNames', ...
T.Properties.VariableNames);
end

hold on;
for i = 1:length(subTables)
    data_per_dist = subTables{i,1};
    x = double(string(data_per_dist.sample_power));
    y = double(data_per_dist.mean_cpu_time);    

    x = log(x)/log(2);
    y = log(y);

    plot(x,y,'.-','DisplayName',string(data_per_dist.name(1)))
end
bp = gca;
bp.XAxis.TickLabelInterpreter = 'latex';
xlabel('$log_{2}(N)$','Interpreter','latex')
ylabel('log(CPU) Time','Interpreter','latex')
legend('Location','northwest')
xlim(xlimits)
ylim(ylimits)
if save_figs
    saveas(bp, fullfile(dir_name, [fig_name, '.png']))
end



