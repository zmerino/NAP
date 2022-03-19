
close all; clear;

filename = 'br0_table.dat';

dirname = fullfile('data','theoretical_threshold');

pathname = fullfile(dirname, filename);

data = readtable(pathname);


disp(data)

num_col = size(data, 2)-1;
num_dist = num_col/2;

dist_index = linspace(2,num_col,num_dist);

names = {'Uniform-Mix', 'Generalized-Pareto', 'Stable',...
    'Tri-Modal-Normal', 'Normal', 'Uniform', 'Beta(0.5,1.5)',...
    'Beta(2,0.5)', 'Beta(0.5,0.5)'}';

labels = {'$2^{8}$', '$2^{9}$', '$2^{10}$', '$2^{11}$', '$2^{12}$'};
label_val = [8,9,10,11,12];

% labels = {'$2^{8}$', '$2^{9}$', '$2^{10}$', '$2^{11}$', '$2^{12}$',...
%     '$2^{13}$', '$2^{14}$', '$2^{15}$', '$2^{16}$', '$2^{17}$',...
%     '$2^{18}$', '$2^{19}$', '$2^{20}$', '$2^{21}$', '$2^{22}$',...
%     '$2^{23}$', '$2^{24}$', '$2^{25}$', '$2^{26}$'};
% label_val = [8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26];

YexpScale = -2;

min_pow = data{1,1};
max_pow = data{end,1};

xT = linspace(2^min_pow,2^max_pow,10000);
T = 2.6741.*xT.^(-0.626) + 0.005;

figure('Name','testBR')
hold on;
cc=lines(num_dist);
x = data{:,1}';
count = 1;
for j = dist_index

    y = data{:,j}';
    stdev = data{:,j+1}';

    curve1 = y + stdev;
    curve2 = y - stdev;

    x2 = [x, fliplr(x)];
    inBetween = [curve1, fliplr(curve2)];
    f = fill(x2, inBetween,cc(count,:));
    set(f,'facealpha',0.25)
    set(f,'edgealpha',0)

    h(count) = plot(x,y,'DisplayName',char(names(count)),'Color',cc(count,:));

    count = count + 1;
end
s = plot(log(xT)/log(2), T, 'k--','DisplayName','$\Gamma$');
bp = gca;
bp.XAxis.TickLabelInterpreter = 'latex';
bp.XTickLabel = labels;
bp.YAxis.Exponent = YexpScale;
xticks(label_val);
% xlabel('$log_{2}(x)$','Interpreter','latex')
xlabel('$x$','Interpreter','latex')
xticklabels(labels)
ylabel('Threshold','Interpreter','latex')
legend([s,h],'Interpreter','latex')