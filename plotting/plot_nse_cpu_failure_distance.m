

addpath("functions/")
addpath("data/")

save_figs = true;

% Import data table
filename = fullfile('data','estimator_meta_data.dat');
filename = 'nse_estimator_meta_data.dat';
data = readtable(filename);

% Define labels for figures
names = {'Uniform-Mix', 'Generalized-Pareto', 'Stable',...
    'Tri-Modal-Normal', 'Normal', 'Uniform', 'Beta(0.5,1.5)',...
    'Beta(2,0.5)', 'Beta(0.5,0.5)'}';

% labels = {'$2^{8}$', '$2^{9}$', '$2^{10}$'};
% label_val = [8,9,10];

labels = {'$2^{8}$', '$2^{9}$', '$2^{10}$', '$2^{11}$', '$2^{12}$',...
    '$2^{13}$', '$2^{14}$', '$2^{15}$'};

label_val = [8,9,10,11,12,13,14,15];
% 
% labels = {'$2^{8}$', '$2^{9}$', '$2^{10}$', '$2^{11}$', '$2^{12}$',...
%     '$2^{13}$', '$2^{14}$', '$2^{15}$', '$2^{16}$', '$2^{17}$', '$2^{18}$', '$2^{19}$', '$2^{20}$'};
% 
% label_val = [8,9,10,11,12,13,14,15,16,17,18,19,20];

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
legend
if save_figs
    saveas(bp, fullfile('figures', [fig_name, '.png']))
end

data_nse = data('NSE'==convertCharsToStrings(data.estimator),:);
data_nmem = data('NMEM'==convertCharsToStrings(data.estimator),:);


fig_name = 'NSE CPU per Distribution';
figure('Name',fig_name)
b = boxchart(log(data_nse.sample_power)/log(2), data_nse.cpu_time, 'GroupByColor',data_nse.distribution);
bp = gca;
bp.XAxis.TickLabelInterpreter = 'latex';
xlabel('$log_{2}(N)$','Interpreter','latex')
ylabel('CPU Time','Interpreter','latex')
legend('Location','northwest')
if save_figs
    saveas(bp, fullfile('figures', [fig_name, '.png']))
end

fig_name = 'NMEM CPU per Distribution';
figure('Name',fig_name)
b = boxchart(log(data_nmem.sample_power)/log(2), data_nmem.cpu_time, 'GroupByColor',data_nmem.distribution);
bp = gca;
bp.XAxis.TickLabelInterpreter = 'latex';
xlabel('$log_{2}(N)$','Interpreter','latex')
ylabel('CPU Time','Interpreter','latex')
legend('Location','northwest')
if save_figs
    saveas(bp, fullfile('figures', [fig_name, '.png']))
end

% -----------------------------------------------
fig_name = 'Max lagrange multiplier';
figure('Name',fig_name)
b = boxchart(log(data.sample_power)/log(2), data.lagrange, 'GroupByColor',data.estimator);
bp = gca;
bp.XAxis.TickLabelInterpreter = 'latex';
xlabel('$log_{2}(N)$','Interpreter','latex')
ylabel('Max Lagrange Multiplier','Interpreter','latex')
legend
if save_figs
    saveas(bp, fullfile('figures', [fig_name, '.png']))
end

fig_name = 'NSE Max lagrange multiplier';
figure('Name',fig_name)
b = boxchart(log(data_nse.sample_power)/log(2), data_nse.lagrange, 'GroupByColor',data_nse.distribution);
bp = gca;
bp.XAxis.TickLabelInterpreter = 'latex';
xlabel('$log_{2}(N)$','Interpreter','latex')
ylabel('Max Lagrange Multiplier','Interpreter','latex')
legend
if save_figs
    saveas(bp, fullfile('figures', [fig_name, '.png']))
end

fig_name = 'NMEM Max lagrange multiplier';
figure('Name',fig_name)
b = boxchart(log(data_nmem.sample_power)/log(2), data_nmem.lagrange, 'GroupByColor',data_nmem.distribution);
bp = gca;
bp.XAxis.TickLabelInterpreter = 'latex';
xlabel('$log_{2}(N)$','Interpreter','latex')
ylabel('Max Lagrange Multiplier','Interpreter','latex')
legend
if save_figs
    saveas(bp, fullfile('figures', [fig_name, '.png']))
end

% -----------------------------------------------
YexpScale = -2.0;

fig_name = 'KL';
figure('Name',fig_name)
b = boxchart(log(data.sample_power)/log(2), data.kl, 'GroupByColor',data.estimator);
bp = gca;
bp.XAxis.TickLabelInterpreter = 'latex';
xlabel('$log_{2}(N)$','Interpreter','latex')
ylabel('KL','Interpreter','latex')
legend
if save_figs
    saveas(bp, fullfile('figures', [fig_name, '.png']))
end

% TODO: Make sure kl/mse data is positive. Why the negative values?
fig_name = 'KL NSE by distribution';
figure('Name',fig_name)
b = boxchart(log(data_nse.sample_power)/log(2), abs(data_nse.kl),...
        'GroupByColor',data_nse.distribution);
% figure('Name','KL NSE by distribution')
% b = boxchart(log(data_nse.sample_power)/log(2), log(data_nse.kl),...
%         'GroupByColor',data_nse.distribution);
bp = gca;
bp.XAxis.TickLabelInterpreter = 'latex';
bp.XTickLabel = labels;
bp.YAxis.Exponent = YexpScale;
xticks(label_val);
xlabel('$log_{2}(N)$','Interpreter','latex')
xticklabels(labels)
ylabel('$KL$','Interpreter','latex')
legend(unique(data_nse.name),'Interpreter','latex')
if save_figs
    saveas(bp, fullfile('figures', [fig_name, '.png']))
end

fig_name = 'KL NMEM by distribution';
figure('Name',fig_name)
b = boxchart(log(data_nmem.sample_power)/log(2), data_nmem.kl,...
        'GroupByColor',data_nmem.distribution);
bp = gca;
bp.XAxis.TickLabelInterpreter = 'latex';
bp.XTickLabel = labels;
bp.YAxis.Exponent = YexpScale;
xticks(label_val);
xlabel('$log_{2}(N)$','Interpreter','latex')
xticklabels(labels)
ylabel('$KL$','Interpreter','latex')
legend(unique(data_nmem.name),'Interpreter','latex')
if save_figs
    saveas(bp, fullfile('figures', [fig_name, '.png']))
end

% TODO: Make sure kl/mse data is positive. Why the negative values?
fig_name = 'LOG KL NSE by distribution';
figure('Name',fig_name)
b = boxchart(log(data_nse.sample_power)/log(2), log(abs(data_nse.kl)),...
        'GroupByColor',data_nse.distribution);
% figure('Name','KL NSE by distribution')
% b = boxchart(log(data_nse.sample_power)/log(2), log(data_nse.kl),...
%         'GroupByColor',data_nse.distribution);
bp = gca;
bp.XAxis.TickLabelInterpreter = 'latex';
bp.XTickLabel = labels;
bp.YAxis.Exponent = YexpScale;
xticks(label_val);
xlabel('$log_{2}(N)$','Interpreter','latex')
xticklabels(labels)
ylabel('$KL$','Interpreter','latex')
legend(unique(data_nse.name),'Interpreter','latex')
if save_figs
    saveas(bp, fullfile('figures', [fig_name, '.png']))
end

fig_name = 'LOG KL NMEM by distribution';
figure('Name',fig_name)
b = boxchart(log(data_nmem.sample_power)/log(2), log(data_nmem.kl),...
        'GroupByColor',data_nmem.distribution);
bp = gca;
bp.XAxis.TickLabelInterpreter = 'latex';
bp.XTickLabel = labels;
bp.YAxis.Exponent = YexpScale;
xticks(label_val);
xlabel('$log_{2}(N)$','Interpreter','latex')
xticklabels(labels)
ylabel('$KL$','Interpreter','latex')
legend(unique(data_nmem.name),'Interpreter','latex')
if save_figs
    saveas(bp, fullfile('figures', [fig_name, '.png']))
end

% -----------------------------------------------
fig_name = 'MSE';
figure('Name',fig_name)
b = boxchart(log(data.sample_power)/log(2), data.mse, 'GroupByColor',data.estimator);
bp = gca;
bp.XAxis.TickLabelInterpreter = 'latex';
xlabel('$log_{2}(N)$','Interpreter','latex')
ylabel('$L_1$','Interpreter','latex')
legend
if save_figs
    saveas(bp, fullfile('figures', [fig_name, '.png']))
end

fig_name = 'MSE NSE by distribution';
figure('Name',fig_name)
b = boxchart(log(data_nse.sample_power)/log(2), data_nse.mse,...
        'GroupByColor',data_nse.distribution);
bp = gca;
bp.XAxis.TickLabelInterpreter = 'latex';
bp.XTickLabel = labels;
bp.YAxis.Exponent = YexpScale;
xticks(label_val);
xlabel('$log_{2}(N)$','Interpreter','latex')
xticklabels(labels)
ylabel('$L_1$','Interpreter','latex')
legend(unique(data_nse.name),'Interpreter','latex')
if save_figs
    saveas(bp, fullfile('figures', [fig_name, '.png']))
end

fig_name = 'MSE NMEM by distribution';
figure('Name',fig_name)
b = boxchart(log(data_nmem.sample_power)/log(2), data_nmem.mse,...
        'GroupByColor',data_nmem.distribution);
bp = gca;
bp.XAxis.TickLabelInterpreter = 'latex';
bp.XTickLabel = labels;
bp.YAxis.Exponent = YexpScale;
xticks(label_val);
xlabel('$log_{2}(N)$','Interpreter','latex')
xticklabels(labels)
ylabel('$L_1$','Interpreter','latex')
legend(unique(data_nmem.name),'Interpreter','latex')
if save_figs
    saveas(bp, fullfile('figures', [fig_name, '.png']))
end

fig_name = 'LOG MSE NSE by distribution';
figure('Name',fig_name)
b = boxchart(log(data_nse.sample_power)/log(2), log(data_nse.mse),...
        'GroupByColor',data_nse.distribution);
bp = gca;
bp.XAxis.TickLabelInterpreter = 'latex';
bp.XTickLabel = labels;
bp.YAxis.Exponent = YexpScale;
xticks(label_val);
xlabel('$log_{2}(N)$','Interpreter','latex')
xticklabels(labels)
ylabel('$L_1$','Interpreter','latex')
legend(unique(data_nse.name),'Interpreter','latex')
if save_figs
    saveas(bp, fullfile('figures', [fig_name, '.png']))
end

fig_name = 'LOG MSE NMEM by distribution';
figure('Name',fig_name)
b = boxchart(log(data_nmem.sample_power)/log(2), log(data_nmem.mse),...
        'GroupByColor',data_nmem.distribution);
bp = gca;
bp.XAxis.TickLabelInterpreter = 'latex';
bp.XTickLabel = labels;
bp.YAxis.Exponent = YexpScale;
xticks(label_val);
xlabel('$log_{2}(N)$','Interpreter','latex')
xticklabels(labels)
ylabel('$L_1$','Interpreter','latex')
legend(unique(data_nmem.name),'Interpreter','latex')
if save_figs
    saveas(bp, fullfile('figures', [fig_name, '.png']))
end








