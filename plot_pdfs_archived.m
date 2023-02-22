clc; clear; close all;

addpath("functions/")

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% To Plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

plot_mse_dist = true;
plot_pdf = true;
plot_cdf = false;
plot_sqr = false;
plot_heavy = false;
save_figs = false;

%%%%%%%%%%%%%%%%%%%%%%%%%%% Import data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Create actual distrobution onjects
actual = distributions;

% Path to directory
dir_name = fullfile('data','pdf_estimates');
dir_name = fullfile('data','cpu_20_t_50_maxN_22_pdf_estimates');
dir_name = fullfile('data','cpu_15_t_50_maxN_22_pdf_estimates');
% cpu_20_t_50_maxN_22_pdf_estimates

% Define the etimates to plot

% Sample range
n_vec = 2.^[8,9,10,11,12,13,14,15,16,17,18];
n_vec = 2.^[14];
% Distribution range
% d_vec = ["Beta-a0p5-b1p5","Beta-a2-b0p5","Beta-a0p5-b0p5","Generalized-Pareto","Stable"];
% d_vec = ["Trimodal-Normal","Uniform","Normal","Uniform-Mix","Beta-a0p5-b1p5","Beta-a2-b0p5","Beta-a0p5-b0p5","Generalized-Pareto","Stable"];
d_vec = ["Trimodal-Normal","Uniform","Normal"];
% d_vec = ["Beta-a0p5-b1p5","Beta-a2-b0p5","Beta-a0p5-b0p5","Generalized-Pareto"];
% d_vec = ["Normal"];

% Trials per sample
trials = 50;

% Initialize empty structure to store data in
n = length(n_vec);
d = length(d_vec);
t = length(d_vec);
nse = cell(d,n,t,7);
nmem = cell(d,n,t,7);

for i = 1:length(d_vec)
    for j = 1:length(n_vec)
        for k = 1:trials
                        
            % Make sure that distribution names are saved as character vectors

            % distribution names
            nse{i,j,k,1} = d_vec(i);
            nmem{i,j,k,1} = d_vec(i);
            % sample size
            nse{i,j,k,2} = n_vec(j);
            nmem{i,j,k,2} = n_vec(j);

            % Load data from estimate directory
            load(fullfile(dir_name,['nse_',num2str(d_vec(i)),'_t_',num2str(trials),'_s_', num2str(n_vec(j)),'.mat']));
            load(fullfile(dir_name,['nmem_',num2str(d_vec(i)),'_t_',num2str(trials),'_s_', num2str(n_vec(j)),'.mat']));

            % x
            nse_data = nse_pdf_data;
            nmem_data = nmem_pdf_data;

            nse{i,j,k,3} = nse_data{2,k}'; % make sure column vector
            nmem{i,j,k,3} = nmem_data{2,k};

            % pdf
            nse{i,j,k,4} = nse_data{3,k}'; % make sure column vector
            nmem{i,j,k,4} = nmem_data{3,k};

            % cdf
            nse{i,j,k,5} = trapz(nse_data{2,k}, nse_data{3,k});
            nmem{i,j,k,5} = trapz(nmem_data{2,k}, nmem_data{3,k});

            % u and sqr
            [nse{i,j,k,6}, nse{i,j,k,7}] = misc_functions.sqr(nse_data{2,k},nse_data{3,k},nse_data{1,k});
            [nmem{i,j,k,6}, nmem{i,j,k,7}] = misc_functions.sqr(nmem_data{2,k},nmem_data{3,k},nmem_data{1,k});



        end
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%% Ploting data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define labels for figures
names = {'Uniform-Mix', 'Generalized-Pareto', 'Stable',...
    'Tri-Modal-Normal', 'Normal', 'Uniform', 'Beta(0.5,1.5)',...
    'Beta(2,0.5)', 'Beta(0.5,0.5)'}';
names = {'Beta(0.5,1.5)','Beta(2,0.5)', 'Beta(0.5,0.5)', 'Generalized-Pareto', 'Stable'}';


%%% MSE Plots %%%
actual = distributions;
actual.generate_data = false;

if plot_mse_dist

    qnum = 20;

    % x values for given quantiles
    for d_idx = 1:size(nse, 1)

        for n_idx = 1:size(nse, 2)

            mse_dist_nse = zeros(qnum,trials);
            mse_dist_nmem = zeros(qnum,trials);
%             for t_idx = 1:size(nse, 3)
            for t_idx = 1:3

                
                disp(['dist: ',num2str(d_idx),'/',num2str(size(nse, 1)),...
                    ' s: ',num2str(n_idx),'/',num2str(size(nse, 2)),...
                    ' t: ',num2str(t_idx),'/',num2str(size(nse, 3))])

                
               % NMEM
                xs = nmem{d_idx,n_idx,t_idx,3}(:,1);
                fs = nmem{d_idx,n_idx,t_idx,4}(:,1);
                Fs = nmem{d_idx,n_idx,t_idx,5}(:,1);

                [metric_dist, quantiles,xq,fq] = quant_metric(xs, fs, Fs, qnum, d_vec(d_idx));
                mse_dist_nmem(:,t_idx) = metric_dist;

                %NSE
                xs = nse{d_idx,n_idx,t_idx,3}(:,1);
                fs = nse{d_idx,n_idx,t_idx,4}(:,1);
                Fs = nse{d_idx,n_idx,t_idx,5}(:,1);

                [metric_dist,~,~,~] = quant_metric(xs, fs, Fs, qnum, d_vec(d_idx));
                mse_dist_nse(:,t_idx) = metric_dist;


%                 q_min = 0.01;
%                 q_max = 0.99;
%                 
%                 xs = nse{d_idx,n_idx,t_idx,3}(:,1);
%                 fs = nse{d_idx,n_idx,t_idx,4}(:,1);
%                 Fs = nse{d_idx,n_idx,t_idx,5}(:,1);


%                 quantiles = linspace(q_min, q_max, qnum);
%                 xq = interp1(Fs,xs,quantiles);
%                 fq = interp1(xs,fs,xq);

%                 % get actual distribution
%                 actual.min_limit = min(xs);
%                 actual.max_limit = max(xs);
%                 actual.dist_name = d_vec(d_idx);
% 
%                 for q = 2:qnum
% 
%                     trunc_xs = linspace(xq(q-1), xq(q), 1000);
%                     trunc_fs = interp1(xs,fs,trunc_xs);
% 
%                     actual.x = linspace(min(trunc_xs),max(trunc_xs),length(trunc_xs));
%                     actual = actual.dist_list();
% 
%                     fact = actual.pdf_y;
% 
%                     mse_dist_nse(q-1,t_idx) = mse(fact, trunc_fs);
% 
%                 end

%                 xs = nmem{d_idx,n_idx,t_idx,3}(:,1);
%                 fs = nmem{d_idx,n_idx,t_idx,4}(:,1);
%                 Fs = nmem{d_idx,n_idx,t_idx,5}(:,1);
% 
% 
%                 quantiles = linspace(q_min, q_max, qnum);
%                 xq = interp1(Fs,xs,quantiles);
%                 fq = interp1(xs,fs,xq);
% 
%                 % get actual distribution
%                 actual.min_limit = min(xs);
%                 actual.max_limit = max(xs);
%                 actual.dist_name = d_vec(d_idx);
% 
%                 for q = 2:qnum
% 
%                     trunc_xs = linspace(xq(q-1), xq(q), 1000);
%                     trunc_fs = interp1(xs,fs,trunc_xs);
% 
%                     actual.x = linspace(min(trunc_xs),max(trunc_xs),length(trunc_xs));
%                     actual = actual.dist_list();
% 
%                     fact = actual.pdf_y;
% 
%                     mse_dist_nmem(q-1,t_idx) = mse(fact, trunc_fs);
% 
%                 end

            end

            % plot quantiles
            figure('Name',['quantile_d_',convertStringsToChars(d_vec(d_idx)),'_s_',num2str(n_vec(n_idx))])
            hold on;
            plot(xs, Fs, '--k')
            plot(xq, quantiles, 'og')
            plot(xq, fq, '.b')
            xlabel('x')
            ylabel('Q(x)')

            figure('Name',['d_',convertStringsToChars(d_vec(d_idx)),'_s_',num2str(n_vec(n_idx))])
            hold on;
            for t = 1:size(mse_dist_nse,2)
                plot(quantiles, mse_dist_nse(:,t), '-r')
                plot(quantiles, mse_dist_nmem(:,t), '-b')
            end
            bp = gca;
            bp.YAxis.Scale ="log";
            xlabel('x')
            ylabel('MSE(x)')


            figure('Name',['avg_d_',convertStringsToChars(d_vec(d_idx)),'_s_',num2str(n_vec(n_idx))])
            hold on;
            plot(quantiles, mean(mse_dist_nse,2), '-r')
            plot(quantiles, mean(mse_dist_nmem,2), '-b')
            xlabel('x')
            ylabel('MSE(x)')

        end
    end



end


%%% PDFs %%%
if plot_pdf
    % Plot data the order of indices corresponds to distribution, sample size,
    % trial. Plot all trials per sample per distrobution
    for d_idx = 1:size(nse, 1)

        for n_idx = 1:size(nse, 2)

            fig_name = ['PDF_d',convertStringsToChars(d_vec(d_idx)),num2str(d_idx),'_s_',num2str(n_vec(n_idx)),'_t_',num2str(trials)];
            figure('Name',fig_name)
            hold on;
            % Set actual distribution type
            actual.dist_name = nse{d_idx,n_idx,1,1}(1,:);
            % Set distrobution over specific range i.e. the x values from the
            % first trial estimate

            actual.x = nse{d_idx,n_idx,1,3}(:,1);
            % Generate distrobution attributes
            actual = dist_list(actual);

            % Remove zeros in actual.pdf_y generated from inf/nan values
            [row,~] = find(actual.pdf_y > 0);
            for s_idx = 1:trials
                nse_h = plot(nse{d_idx,n_idx,s_idx,3}, nse{d_idx,n_idx,s_idx,4},'-r', 'DisplayName', '$\hat{f}_i^{NSE}(x)$');
                nmem_h = plot(nmem{d_idx,n_idx,s_idx,3}, nmem{d_idx,n_idx,s_idx,4},'-b', 'DisplayName', '$\hat{f}_i^{NMEM}(x)$');
            end
            g = plot(actual.x(row,1), actual.pdf_y(row,1), '-k', 'DisplayName','$f(x)$');
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
            if save_figs
                saveas(bp, fullfile('figures', [fig_name, '.png']))
            end
        end
    end
end

%%% CDFs %%%
if plot_cdf
    % Plot data the order of indices corresponds to distribution, sample size,
    % trial. Plot all trials per sample per distrobution
    for d_idx = 1:size(nse, 1)

        for n_idx = 1:size(nse, 2)

            fig_name = ['CDF_d',convertStringsToChars(d_vec(d_idx)),num2str(d_idx),'_s_',num2str(n_vec(n_idx)),'_t_',num2str(trials)];
            figure('Name',fig_name)
            hold on;
            % Set actual distribution type
            actual.dist_name = nse{d_idx,n_idx,1,1};
            % Set distrobution over specific range i.e. the x values from the
            % first trial estimate
            actual.x = nse{d_idx,n_idx,1,3}(:,1);
            % Generate distrobution attributes
            actual = dist_list(actual);

            % Remove zeros in actual.pdf_y generated from inf/nan values
            [row,~] = find(actual.pdf_y > 0);

            for s_idx = 1:trials

                nse_h = plot(nse{d_idx,n_idx,s_idx,3}, nse{d_idx,n_idx,s_idx,5},'-r', 'DisplayName', '$\hat{F}_i^{NSE}(x)$');
                nmem_h = plot(nmem{d_idx,n_idx,s_idx,3}, nmem{d_idx,n_idx,s_idx,5},'-b', 'DisplayName', '$\hat{F}_i^{NMEM}(x)$');
            end
            g = plot(actual.x(row,1), actual.cdf_y(row,1), '-k', 'DisplayName','$F(x)$');
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
            if save_figs
                saveas(bp, fullfile('figures', [fig_name, '.png']))
            end
        end
    end
end

%%% SQRs %%%

if plot_sqr
    % Plot data the order of indices corresponds to distribution, sample size,
    % trial. Plot all trials per sample per distrobution
    for d_idx = 1:size(nse, 1)

        for n_idx = 1:size(nse, 2)

            fig_name = ['SQR_d',convertStringsToChars(d_vec(d_idx)),num2str(d_idx),'_s_',num2str(n_vec(n_idx)),'_t_',num2str(trials)];
            figure('Name',fig_name)
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
            % Set actual distribution type
            actual.dist_name = nse{d_idx,n_idx,1,1};
            % Set distrobution over specific range i.e. the x values from the
            % first trial estimate
            actual.x = nse{d_idx,n_idx,1,3}(:,1);
            % Generate distrobution attributes
            actual = dist_list(actual);

            % Remove zeros in actual.pdf_y generated from inf/nan values
            [row,~] = find(actual.pdf_y > 0);

            for s_idx = 1:trials
                nse_h = plot(nse{d_idx,n_idx,s_idx,6}, nse{d_idx,n_idx,s_idx,7},'-r', 'DisplayName', '$\hat{sqr}_i^{NSE}(x)$');
                nmem_h = plot(nmem{d_idx,n_idx,s_idx,6}, nmem{d_idx,n_idx,s_idx,7},'-b', 'DisplayName', '$\hat{sqr}_i^{NMEM}(x)$');
            end
            bp = gca;
            xlabel('$x$','Interpreter','latex')
            ylabel('$\hat{SQR}(x)$','Interpreter','latex')
            legend([nse_h(1),nmem_h(1)],'Interpreter','latex')
            xlabel('$x$','Interpreter','latex')
            ylabel('$sqr(x)$','Interpreter','latex')
            if save_figs
                saveas(bp, fullfile('figures', [fig_name, '.png']))
            end
        end
    end
end

%%% Heavy tails %%%

if plot_heavy
    % Plot data the order of indices corresponds to distribution, sample size,
    % trial. Plot all trials per sample per distrobution
    for d_idx = 1:size(nse, 1)

        for n_idx = 1:size(nse, 2)

            fig_name = ['Heavy_Tails_d',convertStringsToChars(d_vec(d_idx)),num2str(d_idx),'_s_',num2str(n_vec(n_idx)),'_t_',num2str(trials)];
            figure('Name',fig_name)
            hold on;
            for j = 1:size(nse{d_idx,n_idx,3}, 2)

                % Set actual distribution type
                actual.dist_name = nse{d_idx,n_idx,1,1};
                % Set distrobution over specific range i.e. the x values from the
                % first trial estimate
                xs = nse{d_idx,n_idx,1,3}(:,1);
                x_min = min(xs);
                x_max = max(xs);
                actual.x = linspace(x_min,x_max*2,1000);
                % Generate distrobution attributes
                actual = dist_list(actual);

                % Remove zeros in actual.pdf_y generated from inf/nan values
                [row,~] = find(actual.pdf_y > 0);

                for s_idx = 1:trials
                    [nse_htx, nse_hty, nse_htx_act, nse_hty_act] = heavy(nse{d_idx,n_idx,s_idx,3}(:,j),...
                        nse{d_idx,n_idx,s_idx,5}(:,j), actual);


                    test = nmem{d_idx,n_idx,s_idx,3}(:,j);
                    test2 = nmem{d_idx,n_idx,s_idx,5}(:,j);
                    test3 = nmem{d_idx,n_idx,s_idx,3};
                    test4 = nmem{d_idx,n_idx,s_idx,5};
                    [nmem_htx, nmem_hty, nmem_htx_act, nmem_hty_act] = heavy(nmem{d_idx,n_idx,s_idx,3}(:,j),...
                        nmem{d_idx,n_idx,s_idx,5}(:,j), actual);

                    nse_h = plot(real(nse_htx), real(nse_hty), '-r', 'DisplayName', '$\hat{f}_i^{NSE}(x)$');
                    nmem_h = plot(real(nmem_htx), real(nmem_hty), '-b', 'DisplayName', '$\hat{f}_i^{NMEM}(x)$');
                end
            end
            g = plot(nse_htx_act,nse_hty_act,'-k', 'DisplayName','$f(x)$');
            bp = gca;
            set(bp, 'YScale', 'log')
            xlabel('$x$','Interpreter','latex')
            ylabel('$log(1 - \hat{F}(x))$','Interpreter','latex')
            %             xlim([min(nse_htx),max(nse_htx) + 0.1*(max(nse_htx) - min(nse_htx))])
            %             ylim([min(nse_hty),max(nse_hty) + 0.1*(max(nse_hty) - min(nse_hty))])
            legend([nse_h(1),nmem_h(1), g],'Interpreter','latex', 'Location','southwest')
            if save_figs
                saveas(bp, fullfile('figures', [fig_name, '.png']))
            end
        end
    end
end




