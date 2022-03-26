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

n_vec = 2.^[10,11,12,13,14,15];
% Distribution range
d_vec = ["Stable"];


% Initialize empty structure to store data in
n = length(n_vec);
d = length(d_vec);
nse = cell(d,n,7);
nmem = cell(d,n,7);

for i = 1:length(d_vec)
    for j = 1:length(n_vec)

        % Make sure that distribution names are saved as character vectors

        % distribution names
        nse{i,j,1} = d_vec(i);
        nmem{i,j,1} = d_vec(i);
        % sample size
        nse{i,j,2} = n_vec(j);
        nmem{i,j,2} = n_vec(j);

        % Load data from estimate directory
        load(fullfile(dir_name,['nse_',num2str(d_vec(i)),'_s_', num2str(n_vec(j)),'.mat']));
        load(fullfile(dir_name,['nmem_',num2str(d_vec(i)),'_s_', num2str(n_vec(j)),'.mat']));

        % x
        nse{i,j,3} = nse_data{1,:};
        nmem{i,j,3} = nmem_data{1,:};
        % pdf
        nse{i,j,4} = nse_data{2,:};
        nmem{i,j,4} = nmem_data{2,:};
        % cdf
        nse{i,j,5} = nse_data{3,:};
        nmem{i,j,5} = nmem_data{3,:};
        % u
        nse{i,j,6} = nse_data{4,:};
        nmem{i,j,6} = nmem_data{4,:};
        % sqr
        nse{i,j,7} = nse_data{5,:};
        nmem{i,j,7} = nmem_data{5,:};
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
    for d_indx = 1:size(nse, 1)
        for n_indx = 1:size(nse, 2)
    
            % Set actual distribution type
            actual.dist_name = nse{d_indx,n_indx,1}(1,:);
            % Set distrobution over specific range i.e. the x values from the 
            % first trial estimate

            actual.x = nse{d_indx,n_indx,3}(:,1);
            % Generate distrobution attributes
            actual = dist_list(actual);
    
    
            % Remove zeros in actual.pdf_y generated from inf/nan values
            [row,~] = find(actual.pdf_y > 0);
    
            figure('Name',['Figure_pdf_d',num2str(d_indx),'_s',num2str(n_indx)])
            hold on;
            nse_h = plot(nse{d_indx,n_indx,3}, nse{d_indx,n_indx,4},'-r', 'DisplayName', '$\hat{f}_i^{NSE}(x)$');
            nmem_h = plot(nmem{d_indx,n_indx,3}, nmem{d_indx,n_indx,4},'-b', 'DisplayName', '$\hat{f}_i^{NMEM}(x)$');
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
            legend([nse_h(1),nmem_h(1), g],'Interpreter','latex')
        end
    end
end

%%% CDFs %%%
if plot_cdf
    % Plot data the order of indices corresponds to distribution, sample size,
    % trial. Plot all trials per sample per distrobution
    for d_indx = 1:size(nse, 1)
        for n_indx = 1:size(nse, 2)
    
            % Set actual distribution type
            actual.dist_name = nse{d_indx,n_indx,1};
            % Set distrobution over specific range i.e. the x values from the 
            % first trial estimate
            actual.x = nse{d_indx,n_indx,3}(:,1);
            % Generate distrobution attributes
            actual = dist_list(actual);
    
    
            % Remove zeros in actual.pdf_y generated from inf/nan values
            [row,~] = find(actual.pdf_y > 0);
    
            figure('Name',['Figure_cdf_d',num2str(d_indx),'_s',num2str(n_indx)])
            hold on;
%             h = plot(nse{d_indx,n_indx,3}, nse(d_indx,n_indx).cdf, 'DisplayName', '$\hat{F}_i(x)$');
            nse_h = plot(nse{d_indx,n_indx,3}, nse{d_indx,n_indx,5},'-r', 'DisplayName', '$\hat{F}_i^{NSE}(x)$');
            nmem_h = plot(nmem{d_indx,n_indx,3}, nmem{d_indx,n_indx,5},'-b', 'DisplayName', '$\hat{F}_i^{NMEM}(x)$');
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
            legend([nse_h(1),nmem_h(1), g],'Interpreter','latex')
        end
    end
end

%%% SQRs %%%

if plot_sqr
    % Plot data the order of indices corresponds to distribution, sample size,
    % trial. Plot all trials per sample per distrobution
    for d_indx = 1:size(nse, 1)
        for n_indx = 1:size(nse, 2)
    
            % Set actual distribution type
            actual.dist_name = nse{d_indx,n_indx,1};
            % Set distrobution over specific range i.e. the x values from the 
            % first trial estimate
            actual.x = nse{d_indx,n_indx,3}(:,1);
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
%             h = plot(nse(d_indx,n_indx).u, nse(d_indx,n_indx).sqr, 'DisplayName', '$\hat{sqr}_i(x)$');
            nse_h = plot(nse{d_indx,n_indx,6}, nse{d_indx,n_indx,7},'-r', 'DisplayName', '$\hat{sqr}_i^{NSE}(x)$');
            nmem_h = plot(nmem{d_indx,n_indx,6}, nmem{d_indx,n_indx,7},'-b', 'DisplayName', '$\hat{sqr}_i^{NMEM}(x)$');
            bp = gca;
            xlabel('$x$','Interpreter','latex')
            ylabel('$\hat{SQR}(x)$','Interpreter','latex')
%             legend([nse_h(1),nmem_h(1), g],'Interpreter','latex')
            xlabel('$x$','Interpreter','latex')
            ylabel('$sqr(x)$','Interpreter','latex')
        end
    end
end

%%% Heavy tails %%%

if plot_heavy
    % Plot data the order of indices corresponds to distribution, sample size,
    % trial. Plot all trials per sample per distrobution
    for d_indx = 1:size(nse, 1)
        for n_indx = 1:size(nse, 2)
    
            figure('Name',['Heavy_pdf_d',num2str(d_indx),'_s',num2str(n_indx)])
            hold on;
            for j = 1:size(nse{d_indx,n_indx,3}, 2)  
    
                % Set actual distribution type
                actual.dist_name = nse{d_indx,n_indx,1};
                % Set distrobution over specific range i.e. the x values from the 
                % first trial estimate
                actual.x = nse{d_indx,n_indx,3}(:,1);
                % Generate distrobution attributes
                actual = dist_list(actual);
        
                % Remove zeros in actual.pdf_y generated from inf/nan values
                [row,~] = find(actual.pdf_y > 0);
    
                [nse_htx, nse_hty, nse_htx_act, nse_hty_act] = heavy(nse{d_indx,n_indx,3}(:,j),...
                    nse{d_indx,n_indx,5}(:,j), actual);
                [nmem_htx, nmem_hty, nmem_htx_act, nmem_hty_act] = heavy(nmem{d_indx,n_indx,3}(:,j),...
                    nmem{d_indx,n_indx,5}(:,j), actual);
    
                nse_h = plot(nse_htx, nse_hty, '-r', 'DisplayName', '$\hat{f}_i^{NSE}(x)$');
                nmem_h = plot(nmem_htx, nmem_hty, '-b', 'DisplayName', '$\hat{f}_i^{NMEM}(x)$');

            end
            g = plot(nse_htx_act,nse_hty_act,'--k', 'DisplayName','$f(x)$');
            bp = gca;
            xlabel('$x$','Interpreter','latex')
            ylabel('$\hat{f}(x)$','Interpreter','latex')
            xlim([min(nse_htx),max(nse_htx) + 0.1*(max(nse_htx) - min(nse_htx))])
%             ylim([min(nse_hty),max(nse_hty) + 0.1*(max(nse_hty) - min(nse_hty))])
            legend([nse_h(1),nmem_h(1), g],'Interpreter','latex')
        end
    end
end




