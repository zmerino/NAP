
close all;

addpath("functions/")
addpath("functions_plotting/")
% addpath("data/")
addpath("data_2/")


publicationQuality();


save_figs = true;

% dir_name = fullfile('figures_manuscript','cpu_v1');
dir_name = fullfile('figures_manuscript_final','run_times');
status = mkdir(dir_name);

% Import data table
% filename = fullfile('data_3','meta_data_100.dat');
filename = fullfile('data_test', 'meta_data_4.dat');
filename = fullfile('data_4', 'meta_data_100.dat');

data = readtable(filename);

% Define labels for figures
names = {'Uniform-Mix', 'Generalized-Pareto', 'Stable',...
    'Tri-Modal-Normal', 'Normal', 'Uniform', 'Beta(0.5,1.5)',...
    'Beta(2,0.5)', 'Beta(0.5,0.5)'}';

names = ["Tri-Modal-Normal","Uniform", "Normal", "Beta(0.5,1.5)", "Beta(2,0.5)", "Beta(0.5,0.5)", "Generalized-Pareto"];
% names = ["Tri-Modal-Normal","Uniform", "Normal"];
% names = [ "Normal", "Beta(0.5,0.5)", "Generalized-Pareto"];
% names = ["Beta(0.5,1.5)", "Beta(2,0.5)", "Beta(0.5,0.5)", "Generalized-Pareto"];
names = ["Beta(0.5,1.5)", "Beta(2,0.5)", "Beta(0.5,0.5)"];

labels = {'$2^{8}$', '$2^{9}$', '$2^{10}$', '$2^{11}$', '$2^{12}$',...
    '$2^{13}$', '$2^{14}$', '$2^{15}$', '$2^{16}$', '$2^{17}$', '$2^{18}$'};

label_val = [8,9,10,11,12,13,14,15,16,17,18];

% labels = {'$2^{8}$', '$2^{9}$', '$2^{10}$', '$2^{11}$', '$2^{12}$',...
%     '$2^{13}$', '$2^{14}$', '$2^{15}$', '$2^{16}$', '$2^{17}$', '$2^{18}$', '$2^{19}$', '$2^{20}$'};
% 
% label_val = [8,9,10,11,12,13,14,15,16,17,18,19,20];



names = ["Tri-Modal-Normal","Uniform","Normal","Beta-a0p5-b1p5","Beta-a2-b0p5","Beta-a0p5-b0p5","Generalized-Pareto","Stable"];
names = ["Tri-Modal-Normal","Uniform", "Normal", "Generalized-Pareto"];
labels = {'$2^{8}$', '$2^{9}$', '$2^{10}$', '$2^{11}$', '$2^{12}$',...
    '$2^{13}$', '$2^{14}$', '$2^{15}$', '$2^{16}$', '$2^{17}$',...
    '$2^{18}$', '$2^{19}$', '$2^{20}$', '$2^{21}$', '$2^{22}$'};
label_val = [8,9,10,11,12,13,14,15,16,17,18,19,20,21,22];


% names = [ "Beta(0.5,0.5)"];
% labels = {'$2^{13}$', '$2^{14}$', '$2^{15}$', '$2^{16}$', '$2^{17}$'};
% label_val = [13,14,15,16,17];

% labels = {'$2^{8}$', '$2^{9}$', '$2^{10}$', '$2^{11}$', '$2^{12}$',...
%     '$2^{13}$', '$2^{14}$', '$2^{15}$', '$2^{16}$', '$2^{17}$',...
%     '$2^{18}$', '$2^{19}$', '$2^{20}$', '$2^{21}$', '$2^{22}$',...
%     '$2^{23}$', '$2^{24}$', '$2^{25}$', '$2^{26}$'};
% label_val = [8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26];


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

data_nap_parallel = data('NSE_{parallel}'==convertCharsToStrings(data.estimator),:);
data_nap_serial = data('NSE_{serial}'==convertCharsToStrings(data.estimator),:);
data_nmem = data('NMEM'==convertCharsToStrings(data.estimator),:);

% box plot per distribution
for idx = 1:length(names)
    
    tester = data.name(idx);
    mask = data.name == names(idx);
    data_dist = data(mask,:);
    data_dist.cpu_time = log10(data_dist.cpu_time);
    
    fig_name = sprintf('cpu_%s',names(idx));
    figure('Name',fig_name)
    b = boxchart(log(data_dist.sample_power)/log(2), data_dist.cpu_time, 'GroupByColor',data_dist.estimator);
    bp = gca;
    bp.XAxis.TickLabelInterpreter = 'latex';
    xlabel('$log_{2}(N)$','Interpreter','latex')
    ylabel('CPU Time [$log_{10}(sec)$]','Interpreter','latex')
    legend('Location','northwest')
    title(convertCharsToStrings(names(idx)))
    if save_figs
        saveas(bp, fullfile(dir_name, [fig_name, '.png']))
    end
end


% average line plot per disrubtion
for idx = 1:length(names)
    
    tester = data.name(idx);
    mask = data.name == names(idx);
    data_dist = data(mask,:);
    data_dist.cpu_time = log10(data_dist.cpu_time);
    
    fig_name = sprintf('cpu_%s',names(idx));
    figure('Name',fig_name)
    b = boxchart(log(data_dist.sample_power)/log(2), data_dist.cpu_time, 'GroupByColor',data_dist.estimator);
    bp = gca;
    bp.XAxis.TickLabelInterpreter = 'latex';
    xlabel('$log_{2}(N)$','Interpreter','latex')
    ylabel('CPU Time [$log_{10}(sec)$]','Interpreter','latex')
    legend('Location','northwest')
    title(convertCharsToStrings(names(idx)))
    if save_figs
        saveas(bp, fullfile(dir_name, [fig_name, '.png']))
    end
end

fig_name = 'NAP_{serial} CPU per Distribution';
figure('Name',fig_name)
b = boxchart(log(data_nap_serial.sample_power)/log(2), log(data_nap_serial.cpu_time), 'GroupByColor',data_nap_serial.distribution);
bp = gca;
bp.XAxis.TickLabelInterpreter = 'latex';
xlabel('$log_{2}(N)$','Interpreter','latex')
ylabel('log(CPU) Time','Interpreter','latex')
legend('Location','northwest')
if save_figs
    saveas(bp, fullfile(dir_name, [fig_name, '.png']))
end

fig_name = 'NMEM CPU per Distribution';
figure('Name',fig_name)
b = boxchart(log(data_nmem.sample_power)/log(2), log(data_nmem.cpu_time), 'GroupByColor',data_nmem.distribution);
bp = gca;
bp.XAxis.TickLabelInterpreter = 'latex';
xlabel('$log_{2}(N)$','Interpreter','latex')
ylabel('log(CPU) Time','Interpreter','latex')
legend('Location','northwest')
if save_figs
    saveas(bp, fullfile(dir_name, [fig_name, '.png']))
end

% fig_name = 'NAP_{parallel} Block Size Max';
% figure('Name',fig_name)
% b = boxchart(log(data_nap_parallel.sample_power)/log(2), data_nap_parallel.max_size, 'GroupByColor',data_nap_parallel.distribution);
% bp = gca;
% bp.XAxis.TickLabelInterpreter = 'latex';
% xlabel('$log_{2}(N)$','Interpreter','latex')
% ylabel('$log_{10}(N_b)$','Interpreter','latex')
% legend('Location','northwest')
% bp.YAxis.Scale ="log";
% if save_figs
%     saveas(bp, fullfile(dir_name, [fig_name, '.png']))
% end
% 
% fig_name = 'NAP_{parallel} Block Scale Max';
% figure('Name',fig_name)
% b = boxchart(log(data_nap_parallel.sample_power)/log(2), data_nap_parallel.max_scale, 'GroupByColor',data_nap_parallel.distribution);
% bp = gca;
% bp.XAxis.TickLabelInterpreter = 'latex';
% xlabel('$log_{2}(N)$','Interpreter','latex')
% ylabel('$log_{10}(N_b)$','Interpreter','latex')
% legend('Location','northwest')
% bp.YAxis.Scale ="log";
% if save_figs
%     saveas(bp, fullfile(dir_name, [fig_name, '.png']))
% end


fig_name = 'NAP_{parallel} Block Size Mean';
figure('Name',fig_name)
b = boxchart(log(data_nap_parallel.sample_power)/log(2), data_nap_parallel.mean_size, 'GroupByColor',data_nap_parallel.distribution);
bp = gca;
bp.XAxis.TickLabelInterpreter = 'latex';
xlabel('$log_{2}(N)$','Interpreter','latex')
ylabel('$log_{10}(N_b)$','Interpreter','latex')
legend('Location','northwest')
bp.YAxis.Scale ="log";
if save_figs
    saveas(bp, fullfile(dir_name, [fig_name, '.png']))
end

fig_name = 'NAP_{parallel} Block Scale Mean';
figure('Name',fig_name)
b = boxchart(log(data_nap_parallel.sample_power)/log(2), data_nap_parallel.mean_scale, 'GroupByColor',data_nap_parallel.distribution);
bp = gca;
bp.XAxis.TickLabelInterpreter = 'latex';
xlabel('$log_{2}(N)$','Interpreter','latex')
ylabel('$log_{10}(N_b)$','Interpreter','latex')
legend('Location','northwest')
bp.YAxis.Scale ="log";
if save_figs
    saveas(bp, fullfile(dir_name, [fig_name, '.png']))
end


% -----------------------------------------------
fig_name = 'Max lagrange multiplier';
figure('Name',fig_name)
b = boxchart(log(data.sample_power)/log(2), data.lagrange, 'GroupByColor',data.estimator);
bp = gca;
bp.XAxis.TickLabelInterpreter = 'latex';
xlabel('$log_{2}(N)$','Interpreter','latex')
ylabel('Max Lagrange Multiplier','Interpreter','latex')
legend('Location','northwest')
if save_figs
    saveas(bp, fullfile(dir_name, [fig_name, '.png']))
end


fig_name = 'NAP_{parallel} Avg lagrange multiplier';
figure('Name',fig_name)
b = boxchart(log(data_nap_parallel.sample_power)/log(2), data_nap_parallel.avg_lagrange, 'GroupByColor',data_nap_parallel.distribution);
bp = gca;
bp.XAxis.TickLabelInterpreter = 'latex';
xlabel('$log_{2}(N)$','Interpreter','latex')
ylabel('CPU Time','Interpreter','latex')
legend('Location','northwest')
if save_figs
    saveas(bp, fullfile(dir_name, [fig_name, '.png']))
end

fig_name = 'NMEM Avg lagrange multiplier';
figure('Name',fig_name)
b = boxchart(log(data_nmem.sample_power)/log(2), data_nmem.avg_lagrange, 'GroupByColor',data_nmem.distribution);
bp = gca;
bp.XAxis.TickLabelInterpreter = 'latex';
xlabel('$log_{2}(N)$','Interpreter','latex')
ylabel('CPU Time','Interpreter','latex')
legend('Location','northwest')
if save_figs
    saveas(bp, fullfile(dir_name, [fig_name, '.png']))
end



data_avg_lg = data('NAP_{serial}'~=convertCharsToStrings(data.estimator),:);
fig_name = 'Average vs NMEM lagrange multiplier';
figure('Name',fig_name)
b = boxchart(log(data_avg_lg.sample_power)/log(2), data_avg_lg.avg_lagrange, 'GroupByColor',data_avg_lg.estimator);
bp = gca;
bp.XAxis.TickLabelInterpreter = 'latex';
xlabel('$log_{2}(N)$','Interpreter','latex')
ylabel('Max Lagrange Multiplier','Interpreter','latex')
legend('Location','northwest')
if save_figs
    saveas(bp, fullfile(dir_name, [fig_name, '.png']))
end

fig_name = 'NAP Max lagrange multiplier';
figure('Name',fig_name)
b = boxchart(log(data_nap_parallel.sample_power)/log(2), data_nap_parallel.lagrange, 'GroupByColor',data_nap_parallel.distribution);
bp = gca;
bp.XAxis.TickLabelInterpreter = 'latex';
xlabel('$log_{2}(N)$','Interpreter','latex')
ylabel('Max Lagrange Multiplier','Interpreter','latex')
legend('Location','northwest')
if save_figs
    saveas(bp, fullfile(dir_name, [fig_name, '.png']))
end

fig_name = 'NMEM Max lagrange multiplier';
figure('Name',fig_name)
b = boxchart(log(data_nmem.sample_power)/log(2), data_nmem.lagrange, 'GroupByColor',data_nmem.distribution);
bp = gca;
bp.XAxis.TickLabelInterpreter = 'latex';
xlabel('$log_{2}(N)$','Interpreter','latex')
ylabel('Max Lagrange Multiplier','Interpreter','latex')
legend('Location','northwest')
if save_figs
    saveas(bp, fullfile(dir_name, [fig_name, '.png']))
end



% line plots ==============================================================

xlimits = [7,23];
ylimits = [-4,6];

fig_name = 'line_NAP_{parallel} Wall Time per Distribution';
figure('Name',fig_name)
% compute average ---------------------------------------------------------
data_nap_parallel = convertvars(data_nap_parallel,["distribution","name","estimator","sample_power"],"categorical");
mean_data__parallel = groupsummary(data_nap_parallel,["distribution","name","estimator","sample_power"],"mean",["wall_time","cpu_time"]);
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
    data = subTables{i};
    x = double(string(data.sample_power));
    y = double(data.mean_wall_time);    

    x = log(x)/log(2);
    y = log(y);

    plot(x,y,'.-','DisplayName',string(data.name(1)))
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

fig_name = 'line_NMEM Wall Time per Distribution';
figure('Name',fig_name)
% compute average ---------------------------------------------------------
data_nmem = convertvars(data_nmem,["distribution","name","estimator","sample_power"],"categorical");
mean_data_nmem = groupsummary(data_nmem,["distribution","name","estimator","sample_power"],"mean",["wall_time","cpu_time"]);
T = mean_data_nmem;
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
    data = subTables{i};
    x = double(string(data.sample_power));
    y = double(data.mean_wall_time);    

    x = log(x)/log(2);
    y = log(y);

    plot(x,y,'.-', 'DisplayName',string(data.name(1)))

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


xlimits = [7,23];
ylimits = [-4,8];

fig_name = 'line_NAP_{parallel} CPU per Distribution';
figure('Name',fig_name)
% compute average ---------------------------------------------------------
data_nap_parallel = convertvars(data_nap_parallel,["distribution","name","estimator","sample_power"],"categorical");
mean_data__parallel = groupsummary(data_nap_parallel,["distribution","name","estimator","sample_power"],"mean",["wall_time","cpu_time"]);
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
    data = subTables{i};
    x = double(string(data.sample_power));
    y = double(data.mean_cpu_time);    

    x = log(x)/log(2);
    y = log(y);

    plot(x,y,'.-','DisplayName',string(data.name(1)))
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

fig_name = 'line_NMEM CPU per Distribution';
figure('Name',fig_name)
% compute average ---------------------------------------------------------
data_nmem = convertvars(data_nmem,["distribution","name","estimator","sample_power"],"categorical");
mean_data_nmem = groupsummary(data_nmem,["distribution","name","estimator","sample_power"],"mean",["wall_time","cpu_time"]);
T = mean_data_nmem;
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
    data = subTables{i};
    x = double(string(data.sample_power));
    y = double(data.mean_cpu_time);    

    x = log(x)/log(2);
    y = log(y);

    plot(x,y,'.-', 'DisplayName',string(data.name(1)))

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




% subplots ----------------------------------------------------------------


