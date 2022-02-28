clc; clear; close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%% Import data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Path to directory
dir_name = fullfile('data','estimates');

% Define the etimates to plot

% Sample range
n_vec = [256,512];
% Distribution range
d_vec = ["Beta-a0p5-b1p5","Beta-a2-b0p5","Beta-a0p5-b0p5",...
        "Stable","Generalized-Pareto"];

% Initialize empty structure to store data in
n = length(n_vec);
d = length(d_vec);
nse_data(d,n) = struct( 'dist', [],...
                        'sample', [],...
                        'x', [],...
                        'pdf', [],...
                        'cdf', [],...
                        'u', [],...
                        'sqr', []);

% Read in data
for i = 1:length(d_vec)
    for j = 1:length(n_vec)

        % Make sure that distribution names are saved as character vectors
        nse_data(i,j).dist = convertStringsToChars(d_vec(i));
        nse_data(i,j).sample = n_vec(j);

        nse_data(i,j).x = readmatrix(fullfile(dir_name,...
            ['nse_x_',num2str(d_vec(i)),'_s_', num2str(n_vec(j)),'.txt']));
        nse_data(i,j).pdf = readmatrix(fullfile(dir_name,...
            ['nse_pdf_',num2str(d_vec(i)),'_s_', num2str(n_vec(j)),'.txt']));
        nse_data(i,j).cdf = readmatrix(fullfile(dir_name,...
            ['nse_cdf_',num2str(d_vec(i)),'_s_', num2str(n_vec(j)),'.txt']));

        nse_data(i,j).u = readmatrix(fullfile(dir_name,...
            ['nse_u_',num2str(d_vec(i)),'_s_', num2str(n_vec(j)),'.txt']));
        nse_data(i,j).sqr = readmatrix(fullfile(dir_name,...
            ['nse_sqr_',num2str(d_vec(i)),'_s_', num2str(n_vec(j)),'.txt']));
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%% Ploting data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define labels for figures
names = {'Uniform-Mix', 'Generalized-Pareto', 'Stable',...
    'Tri-Modal-Normal', 'Normal', 'Uniform', 'Beta(0.5,1.5)',...
    'Beta(2,0.5)', 'Beta(0.5,0.5)'}';


%%% PDFs %%%

% Plot data the order of indices corresponds to distribution, sample size,
% trial
figure('Name','pdf')
hold on;
plot(nse_data(1,1).x(:,1), nse_data(1,1).pdf(:,1),...
    'DisplayName',...
    [nse_data(1,1).dist, '\_s\_',num2str(nse_data(1,1).sample)])
plot(nse_data(1,1).x(:,2), nse_data(1,1).pdf(:,2),...
    'DisplayName',...
    [nse_data(1,1).dist, '\_s\_',num2str(nse_data(1,1).sample)])
plot(nse_data(1,2).x(:,1), nse_data(1,2).pdf(:,1),...
    'DisplayName',...
    [nse_data(1,2).dist,'\_s\_',num2str(nse_data(1,2).sample)])
plot(nse_data(1,2).x(:,2), nse_data(1,2).pdf(:,2),...
    'DisplayName',...
    [nse_data(1,2).dist,'\_s\_',num2str(nse_data(1,2).sample)])
bp = gca;
xlabel('$x$','Interpreter','latex')
ylabel('$\hat{f}(x)$','Interpreter','latex')


%%% CDFs %%%

% Plot data the order of indices corresponds to distribution, sample size,
% trial
figure('Name','cdf')
hold on;
plot(nse_data(1,1).x(:,1), nse_data(1,1).cdf(:,1),...
    'DisplayName',...
    [nse_data(1,1).dist, '\_s\_',num2str(nse_data(1,1).sample)])
plot(nse_data(1,1).x(:,2), nse_data(1,1).cdf(:,2),...
    'DisplayName',...
    [nse_data(1,1).dist, '\_s\_',num2str(nse_data(1,1).sample)])
plot(nse_data(1,2).x(:,1), nse_data(1,2).cdf(:,1),...
    'DisplayName',...
    [nse_data(1,2).dist,'\_s\_',num2str(nse_data(1,2).sample)])
plot(nse_data(1,2).x(:,2), nse_data(1,2).cdf(:,2),...
    'DisplayName',...
    [nse_data(1,2).dist,'\_s\_',num2str(nse_data(1,2).sample)])
bp = gca;
xlabel('$x$','Interpreter','latex')
ylabel('$\hat{F}(x)$','Interpreter','latex')


%%% SQRs %%%

% Plot data the order of indices corresponds to distribution, sample size,
% trial
figure('Name','sqr')
hold on;
smallN = 256;
smallN2 = 258;
graymax = 240;
range = 0:1/(smallN+1):1;
muLD = range*(smallN + 1) / (smallN + 1);
lemonDrop = sqrt(muLD.*(1-muLD)) * 3.4;
sampleCount2 = (smallN + 2):-1:1;
colorRange = (255-graymax)*sampleCount2/(smallN + 2);
base = repmat(graymax, smallN + 2, 1);
col = (base + colorRange') / 255;
rgb = [col col col];
count2 = 1;
for ii = ceil(smallN2/2):smallN2-1
    ix = [ii ii+1 smallN2-ii smallN2-ii+1];
    fill(range(ix), lemonDrop(ix), rgb(count2, :),'edgecolor','none')
    fill(range(ix), -lemonDrop(ix), rgb(count2, :),'edgecolor','none')
    count2 = count2 + 2;
end
plot(muLD,lemonDrop,'k--');
plot(muLD,-lemonDrop,'k--');
plot(nse_data(1,1).u(:,1), nse_data(1,1).sqr(:,1),...
    'DisplayName',...
    [nse_data(1,1).dist, '\_s\_',num2str(nse_data(1,1).sample)])
plot(nse_data(1,1).u(:,2), nse_data(1,1).sqr(:,2),...
    'DisplayName',...
    [nse_data(1,1).dist, '\_s\_',num2str(nse_data(1,1).sample)])
plot(nse_data(1,2).u(:,1), nse_data(1,2).sqr(:,1),...
    'DisplayName',...
    [nse_data(1,2).dist,'\_s\_',num2str(nse_data(1,2).sample)])
plot(nse_data(1,2).u(:,2), nse_data(1,2).sqr(:,2),...
    'DisplayName',...
    [nse_data(1,2).dist,'\_s\_',num2str(nse_data(1,2).sample)])
bp = gca;
xlabel('$x$','Interpreter','latex')
ylabel('$sqr(x)$','Interpreter','latex')

%%% Heavy tails %%%

actual = distributions;

temp_min_limit =            0; %<---- set upper limit for both
actual.min_limit =          temp_min_limit;  %<--- lower limit to plot
temp_max_limit =            1000; %<---- set upper limit for both
actual.max_limit =          temp_max_limit; %<--- upper limit to plot

% find any of the strings in "str" inside of "distribtuionVector"
flag = isbeta(nse_data(1,1).dist);

% Define plot vector for dist_list from 0-1
if flag(j)
    actual.min_limit = 0;
    actual.max_limit = 1;
else
    % Define plot vector for distribution from actual.min_limit-actual.max_limit
    actual.min_limit = 0;
    actual.max_limit = 10;
end
% Current distribution name
actual.dist_name = nse_data(1,1).dist;
% file name for actual distribution. "A_" puts at the top of the folder.
actual.filename = sprintf(['A_', char(actual.dist_name),'_Act']);




[htx, hty, htx_act, hty_act] = heavy(nse_data(1,1).x(:,1),...
    nse_data(1,1).cdf(:,1), actual);

figure('Name','Heavy Tails')
hold on;
plot(htx,hty,'-r')
plot(htx_act,hty_act,'--k')




