
% Import data table
filename = fullfile('data','estimator_meta_data.dat');
data = readtable(filename);

% Define labels for figures
names = {'Uniform-Mix', 'Generalized-Pareto', 'Stable',...
    'Tri-Modal-Normal', 'Normal', 'Uniform', 'Beta(0.5,1.5)',...
    'Beta(2,0.5)', 'Beta(0.5,0.5)'}';

labels = {'$2^{8}$', '$2^{9}$', '$2^{10}$', '$2^{11}$', '$2^{12}$'};
label_val = [8,9,10,11,12];

labels = {'$2^{8}$', '$2^{9}$', '$2^{10}$', '$2^{11}$', '$2^{12}$',...
    '$2^{13}$', '$2^{14}$', '$2^{15}$'};

label_val = [8,9,10,11,12,13,14,15];

% labels = {'$2^{8}$', '$2^{9}$', '$2^{10}$', '$2^{11}$', '$2^{12}$',...
%     '$2^{13}$', '$2^{14}$', '$2^{15}$', '$2^{16}$', '$2^{17}$',...
%     '$2^{18}$', '$2^{19}$', '$2^{20}$', '$2^{21}$', '$2^{22}$',...
%     '$2^{23}$', '$2^{24}$', '$2^{25}$', '$2^{26}$'};
% label_val = [8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26];



figure('Name','CPU Time')
b = boxchart(log(data.sample_power)/log(2), data.cpu_time, 'GroupByColor',data.estimator);
bp = gca;
bp.XAxis.TickLabelInterpreter = 'latex';
xlabel('$log_{2}(N)$','Interpreter','latex')
ylabel('CPU Time','Interpreter','latex')
legend

figure('Name','Failure Rate')
b = boxchart(log(data.sample_power)/log(2), data.fail, 'GroupByColor',data.estimator);
bp = gca;
bp.XAxis.TickLabelInterpreter = 'latex';
xlabel('$log_{2}(N)$','Interpreter','latex')
ylabel('Failure Rate','Interpreter','latex')
legend

data_nse = data('NSE'==convertCharsToStrings(data.estimator),:);
data_nmem = data('NMEM'==convertCharsToStrings(data.estimator),:);

% -----------------------------------------------
figure('Name','Max lagrange multiplier')
b = boxchart(log(data.sample_power)/log(2), data.lagrange, 'GroupByColor',data.estimator);
bp = gca;
bp.XAxis.TickLabelInterpreter = 'latex';
xlabel('$log_{2}(N)$','Interpreter','latex')
ylabel('Max Lagrange Multiplier','Interpreter','latex')
legend

figure('Name','NSE Max lagrange multiplier')
b = boxchart(log(data_nse.sample_power)/log(2), data_nse.lagrange, 'GroupByColor',data_nse.distribution);
bp = gca;
bp.XAxis.TickLabelInterpreter = 'latex';
xlabel('$log_{2}(N)$','Interpreter','latex')
ylabel('Max Lagrange Multiplier','Interpreter','latex')
legend

figure('Name','NMEM Max lagrange multiplier')
b = boxchart(log(data_nmem.sample_power)/log(2), data_nmem.lagrange, 'GroupByColor',data_nmem.distribution);
bp = gca;
bp.XAxis.TickLabelInterpreter = 'latex';
xlabel('$log_{2}(N)$','Interpreter','latex')
ylabel('Max Lagrange Multiplier','Interpreter','latex')
legend

% -----------------------------------------------
YexpScale = -2.0;

figure('Name','KL')
b = boxchart(log(data.sample_power)/log(2), data.kl, 'GroupByColor',data.estimator);
bp = gca;
bp.XAxis.TickLabelInterpreter = 'latex';
xlabel('$log_{2}(N)$','Interpreter','latex')
ylabel('KL','Interpreter','latex')
legend

% TODO: Make sure kl/mse data is positive. Why the negative values?
figure('Name','KL NSE by distribution')
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

figure('Name','KL NMEM by distribution')
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

% -----------------------------------------------
figure('Name','MSE')
b = boxchart(log(data.sample_power)/log(2), data.mse, 'GroupByColor',data.estimator);
bp = gca;
bp.XAxis.TickLabelInterpreter = 'latex';
xlabel('$log_{2}(N)$','Interpreter','latex')
ylabel('$L_1$','Interpreter','latex')
legend

figure('Name','MSE NSE by distribution')
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

figure('Name','MSE NMEM by distribution')
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









