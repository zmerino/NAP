clc; clear; close all;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%% To Plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

plot_pdf = true;
plot_cdf = true;
plot_sqr = true;
plot_heavy = true;

%%%%%%%%%%%%%%%%%%%%%%%%%%% Import data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Create actual distrobution onject
actual = distributions;

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
if plot_pdf
    % Plot data the order of indices corresponds to distribution, sample size,
    % trial. Plot all trials per sample per distrobution
    for d_indx = 1:size(nse_data, 1)
        for n_indx = 1:size(nse_data, 2)
    
            % Set actual distribution type
            actual.dist_name = nse_data(d_indx,n_indx).dist;
            % Set distrobution over specific range i.e. the x values from the 
            % first trial estimate
            actual.x = nse_data(d_indx,n_indx).x(:,1);
            % Generate distrobution attributes
            actual = dist_list(actual);
    
    
            % Remove zeros in actual.pdf_y generated from inf/nan values
            [row,~] = find(actual.pdf_y > 0);
    
            figure('Name',['Figure_pdf_d',num2str(d_indx),'_s',num2str(n_indx)])
            hold on;
            h = plot(nse_data(d_indx,n_indx).x, nse_data(d_indx,n_indx).pdf, 'DisplayName', '$\hat{f}_i(x)$');
            g = plot(actual.x(row,1), actual.pdf_y(row,1), '--k', 'DisplayName','$f(x)$');
            bp = gca;
            xlabel('$x$','Interpreter','latex')
            ylabel('$\hat{f}(x)$','Interpreter','latex')
            if ismember(actual.dist_name, ["Beta-a0p5-b1p5","Beta-a2-b0p5","Beta-a0p5-b0p5"])
                xlim([0,1])
                ylim([0,10])
            else
                xlim([0,10])
                ylim([0,1])
            end
            legend([h(1), g],'Interpreter','latex')
        end
    end
end

%%% CDFs %%%
if plot_cdf
    % Plot data the order of indices corresponds to distribution, sample size,
    % trial. Plot all trials per sample per distrobution
    for d_indx = 1:size(nse_data, 1)
        for n_indx = 1:size(nse_data, 2)
    
            % Set actual distribution type
            actual.dist_name = nse_data(d_indx,n_indx).dist;
            % Set distrobution over specific range i.e. the x values from the 
            % first trial estimate
            actual.x = nse_data(d_indx,n_indx).x(:,1);
            % Generate distrobution attributes
            actual = dist_list(actual);
    
    
            % Remove zeros in actual.pdf_y generated from inf/nan values
            [row,~] = find(actual.pdf_y > 0);
    
            figure('Name',['Figure_cdf_d',num2str(d_indx),'_s',num2str(n_indx)])
            hold on;
            h = plot(nse_data(d_indx,n_indx).x, nse_data(d_indx,n_indx).cdf, 'DisplayName', '$\hat{F}_i(x)$');
            g = plot(actual.x(row,1), actual.cdf_y(row,1), '--k', 'DisplayName','$F(x)$');
            bp = gca;
            xlabel('$x$','Interpreter','latex')
            ylabel('$\hat{F}(x)$','Interpreter','latex')
            if ismember(actual.dist_name, ["Beta-a0p5-b1p5","Beta-a2-b0p5","Beta-a0p5-b0p5"])
                xlim([0,1])
            else
                xlim([0,10])
            end
            ylim([0,1])
            legend([h(1), g],'Interpreter','latex')
        end
    end
end

%%% SQRs %%%

if plot_sqr
    % Plot data the order of indices corresponds to distribution, sample size,
    % trial. Plot all trials per sample per distrobution
    for d_indx = 1:size(nse_data, 1)
        for n_indx = 1:size(nse_data, 2)
    
            % Set actual distribution type
            actual.dist_name = nse_data(d_indx,n_indx).dist;
            % Set distrobution over specific range i.e. the x values from the 
            % first trial estimate
            actual.x = nse_data(d_indx,n_indx).x(:,1);
            % Generate distrobution attributes
            actual = dist_list(actual);
    
    
            % Remove zeros in actual.pdf_y generated from inf/nan values
            [row,~] = find(actual.pdf_y > 0);
    
            figure('Name',['Figure_sqr_d',num2str(d_indx),'_s',num2str(n_indx)])
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
            h = plot(nse_data(d_indx,n_indx).u, nse_data(d_indx,n_indx).sqr, 'DisplayName', '$\hat{sqr}_i(x)$');
            bp = gca;
            xlabel('$x$','Interpreter','latex')
            ylabel('$\hat{SQR}(x)$','Interpreter','latex')
            legend([h(1), g],'Interpreter','latex')
            xlabel('$x$','Interpreter','latex')
            ylabel('$sqr(x)$','Interpreter','latex')
        end
    end
end

%%% Heavy tails %%%

if plot_heavy
    % Plot data the order of indices corresponds to distribution, sample size,
    % trial. Plot all trials per sample per distrobution
    for d_indx = 1:size(nse_data, 1)
        for n_indx = 1:size(nse_data, 2)
    
            figure('Name',['Heavy_pdf_d',num2str(d_indx),'_s',num2str(n_indx)])
            hold on;
            for j = 1:size(nse_data(d_indx,n_indx).x, 2)  
    
                % Set actual distribution type
                actual.dist_name = nse_data(d_indx,n_indx).dist;
                % Set distrobution over specific range i.e. the x values from the 
                % first trial estimate
                actual.x = nse_data(d_indx,n_indx).x(:,1);
                % Generate distrobution attributes
                actual = dist_list(actual);
        
                % Remove zeros in actual.pdf_y generated from inf/nan values
                [row,~] = find(actual.pdf_y > 0);
    
                [htx, hty, htx_act, hty_act] = heavy(nse_data(d_indx,n_indx).x(:,j),...
                    nse_data(d_indx,n_indx).cdf(:,j), actual);
    
                h = plot(htx, hty, '-r', 'DisplayName', '$\hat{f}_i(x)$');
            end
            g = plot(htx_act,hty_act,'--k', 'DisplayName','$f(x)$');
            bp = gca;
            xlabel('$x$','Interpreter','latex')
            ylabel('$\hat{f}(x)$','Interpreter','latex')
            xlim([min(htx),max(htx) + 0.1*(max(htx) - min(htx))])
            ylim([min(hty),max(hty) + 0.1*(max(hty) - min(hty))])
            legend([h(1), g],'Interpreter','latex')
        end
    end
end




