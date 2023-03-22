
close all;

addpath("functions/")
% addpath("data/")
addpath("data_2/")

save_figs = true;

dir_name = fullfile('figures_manuscript','cpu_v1');
status = mkdir(dir_name);

% Import data table
filename = fullfile('data_3','meta_data_100.dat');

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
names = ["Tri-Modal-Normal","Uniform", "Normal", "Beta(0.5,1.5)", "Beta(2,0.5)", "Beta(0.5,0.5)", "Generalized-Pareto","Stable"];
labels = {'$2^{8}$', '$2^{9}$', '$2^{10}$', '$2^{11}$', '$2^{12}$',...
    '$2^{13}$', '$2^{14}$', '$2^{15}$', '$2^{16}$', '$2^{17}$',...
    '$2^{18}$', '$2^{19}$', '$2^{20}$', '$2^{21}$', '$2^{22}$'};
label_val = [8,9,10,11,12,13,14,15,16,17,18,19,20,21,22];

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
%     bp.YAxis.Scale ="log";
    if save_figs
        saveas(bp, fullfile(dir_name, [fig_name, '.png']))
    end
end



fig_name = 'NAP_{parallel} CPU per Distribution';
figure('Name',fig_name)
b = boxchart(log(data_nap_parallel.sample_power)/log(2), log(data_nap_parallel.cpu_time), 'GroupByColor',data_nap_parallel.distribution);
bp = gca;
bp.XAxis.TickLabelInterpreter = 'latex';
xlabel('$log_{2}(N)$','Interpreter','latex')
ylabel('log(CPU) Time','Interpreter','latex')
legend('Location','northwest')
if save_figs
    saveas(bp, fullfile(dir_name, [fig_name, '.png']))
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









