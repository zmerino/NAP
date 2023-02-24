clc; clear; close all;

addpath("functions/")

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% To Plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

plot_mse_dist = true;
table_name = 'mse_kl_set0.dat';

plot_pdf = true;
plot_cdf = false;
plot_sqr = true;
plot_heavy = true;
save_figs = true;

% figure directory
% fig_dir = fullfile('figures','mse_kl_set0');
fig_dir = fullfile('figures','mse_kl_set_qstep_10');
status = mkdir(fig_dir);

% choose to visualise figures or not
fig_plot = 'on';

%%%%%%%%%%%%%%%%%%%%%%%%%%% Import data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Create actual distrobution onjects
actual = distributions;

% Path to directory
% dir_name = fullfile('data','pdf_estimates');
dir_name = fullfile('data');

% Define the etimates to plot

% Sample range
% n_vec = 2.^[8,9,10,11,12,13,14,15,16,17,18];
n_vec = 2.^[14];
n_vec = 2.^[8,9,10,11,12,13,14,15,16,17,18,19,20,21,22];
n_vec = 2.^[14,15];
% Distribution range
% d_vec = ["Beta-a0p5-b1p5","Beta-a2-b0p5","Beta-a0p5-b0p5","Generalized-Pareto","Stable"];
% d_vec = ["Trimodal-Normal","Uniform","Normal","Uniform-Mix","Beta-a0p5-b1p5","Beta-a2-b0p5","Beta-a0p5-b0p5","Generalized-Pareto","Stable"];
d_vec = ["Trimodal-Normal","Uniform","Normal"];
% d_vec = ["Beta-a0p5-b1p5","Beta-a2-b0p5","Beta-a0p5-b0p5","Generalized-Pareto"];
% d_vec = ["Normal"];

% Define labels for figures
% names = {'Uniform-Mix', 'Generalized-Pareto', 'Stable',...
%     'Tri-Modal-Normal', 'Normal', 'Uniform', 'Beta(0.5,1.5)',...
%     'Beta(2,0.5)', 'Beta(0.5,0.5)'}';
% names = {'Beta(0.5,1.5)','Beta(2,0.5)', 'Beta(0.5,0.5)', 'Generalized-Pareto', 'Stable'}';
names = ["Trimodal-Normal","Uniform","Normal"];

% Trials per sample
trials = 50;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Load data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

filename = 'mse_kl_set0.mat';

load(fullfile(dir_name,['nse_',filename]));
load(fullfile(dir_name,['nmem_',filename]));


%%%%%%%%%%%%%%%%%%%%%%%%%%% Ploting data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% MSE Plots %%%
actual = distributions;
actual.generate_data = false;
% tracks kl/mse for entire sample
global_table = table();

if plot_mse_dist
    
    plt_trials = 10;
    qnum = 10;

    % x values for given quantiles
    for d_idx = 1:size(nse, 1)

        % define distribution name for plotting and building a data table
        actual.dist_name = d_vec(d_idx);

        mse_dist_nmem = zeros(size(nse, 2),plt_trials);
        kl_dist_nmem = zeros(size(nse, 2),plt_trials);
        mse_dist_nse = zeros(size(nse, 2),plt_trials);
        kl_dist_nse = zeros(size(nse, 2),plt_trials);
        for n_idx = 1:size(nse, 2)

            mse_dist_nse_per_q = zeros(qnum,plt_trials);
            mse_dist_nmem_per_q = zeros(qnum,plt_trials);
            kl_dist_nmem_per_q = zeros(qnum,plt_trials);
            kl_dist_nse_per_q = zeros(qnum,plt_trials);
            for t_idx = 1:plt_trials

                % track calculations
                disp(['dist: ',num2str(d_idx),'/',num2str(size(nse, 1)),...
                    ' s: ',num2str(n_idx),'/',num2str(size(nse, 2)),...
                    ' t: ',num2str(t_idx),'/',num2str(plt_trials)])

               % NMEM ------------------------------------------------------
                xs = nmem{d_idx,n_idx,t_idx,3}(:,1);
                fs = nmem{d_idx,n_idx,t_idx,4}(:,1);
                Fs = nmem{d_idx,n_idx,t_idx,5}(:,1);
                % get actual distribution
                actual.min_limit = min(xs);
                actual.max_limit = max(xs);
                actual.x = xs;
                actual = actual.dist_list();

               % MSE per quantile
                [mse_dist_nmem_per_q(:,t_idx), quantiles,xq,fq] = quant_metric(xs, fs, Fs, qnum, d_vec(d_idx),'MSE');
                % MSE per trial
                mse_dist_nmem(n_idx, t_idx) = mse(actual.pdf_y, fs);
                % KL per quantile
                [kl_dist_nmem_per_q(:,t_idx), ~,~,~] = quant_metric(xs, fs, Fs, qnum, d_vec(d_idx),'KL');
                % KL per trial
                kl_dist_nmem(n_idx, t_idx) = KLDiv(actual.pdf_y, fs');

               % NSE ------------------------------------------------------

                xs = nse{d_idx,n_idx,t_idx,3}(:,1);
                fs = nse{d_idx,n_idx,t_idx,4}(:,1);
                Fs = nse{d_idx,n_idx,t_idx,5}(:,1);
                % get actual distribution
                actual.min_limit = min(xs);
                actual.max_limit = max(xs);
                actual.x = xs;
                actual = actual.dist_list();

               % MSE per quantile
                [mse_dist_nse_per_q(:,t_idx),~,~,~] = quant_metric(xs, fs, Fs, qnum, d_vec(d_idx),'MSE');
                % MSE per trial
                mse_dist_nse(n_idx, t_idx) = mse(actual.pdf_y, fs);
                % KL per quantile
                [kl_dist_nse_per_q(:,t_idx), ~,~,~] = quant_metric(xs, fs, Fs, qnum, d_vec(d_idx),'KL');
                % KL per trial
                kl_dist_nse(n_idx, t_idx) = KLDiv(actual.pdf_y, fs');


            end

            % plot quantiles
            fig_name = ['quantile_d_',convertStringsToChars(d_vec(d_idx)),'_s_',num2str(n_vec(n_idx))];
            figure('Name',fig_name, 'visible',fig_plot)
            hold on;
            act_h = plot(xs, Fs, '-k', 'DisplayName', '$F(x)$');
            q_h = plot(xq, quantiles, 'og', 'DisplayName', sprintf('$Q_{%s}(x)$',[num2str(100/qnum),'\%']));
            pdf_h = plot(xq, fq, '.b', 'DisplayName', '$f(Q^{-1}(y))$');
            xlabel('x')
            ylabel('Q(x)')
            legend([act_h(1),q_h(1),pdf_h(1)],'Interpreter','latex', 'Location','northwest')
            bp = gca;
            if save_figs
                saveas(bp, fullfile(fig_dir, [fig_name, '.png']))
            end

            fig_name = ['mse_d_',convertStringsToChars(d_vec(d_idx)),'_s_',num2str(n_vec(n_idx))];
            figure('Name',fig_name, 'visible',fig_plot)
            hold on;
            for t = 1:size(mse_dist_nse_per_q,2)
                nse_h = plot(quantiles, mse_dist_nse_per_q(:,t), '-r', 'DisplayName', '$\hat{f}_i^{NSE}(x)$');
                nmem_h = plot(quantiles, mse_dist_nmem_per_q(:,t), '-b', 'DisplayName', '$\hat{f}_i^{NMEM}(x)$');
            end
            bp = gca;
%             bp.YAxis.Scale ="log";
            xlabel('x')
            ylabel('MSE(x)')
            legend([nse_h(1),nmem_h(1)],'Interpreter','latex')
            if save_figs
                saveas(bp, fullfile(fig_dir, [fig_name, '.png']))
            end


            fig_name = ['avg_mse_d_',convertStringsToChars(d_vec(d_idx)),'_s_',num2str(n_vec(n_idx))];
            figure('Name',fig_name, 'visible',fig_plot)
            hold on;
            nse_h = plot(quantiles, mean(mse_dist_nse_per_q,2), '-r', 'DisplayName', '$\hat{f}_i^{NSE}(x)$');
            nmem_h = plot(quantiles, mean(mse_dist_nmem_per_q,2), '-b', 'DisplayName', '$\hat{f}_i^{NMEM}(x)$');
            xlabel('x')
            ylabel('MSE(x)')
            legend([nse_h(1),nmem_h(1)],'Interpreter','latex')
            bp = gca;
            if save_figs
                saveas(bp, fullfile(fig_dir, [fig_name, '.png']))
            end

            fig_name = ['kl_d_',convertStringsToChars(d_vec(d_idx)),'_s_',num2str(n_vec(n_idx))];
            figure('Name',fig_name, 'visible',fig_plot)
            hold on;
            for t = 1:size(kl_dist_nmem_per_q,2)
                nse_h = plot(quantiles, kl_dist_nse_per_q(:,t), '-r', 'DisplayName', '$\hat{f}_i^{NSE}(x)$');
                nmem_h = plot(quantiles, kl_dist_nmem_per_q(:,t), '-b', 'DisplayName', '$\hat{f}_i^{NMEM}(x)$');
            end
            bp = gca;
%             bp.YAxis.Scale ="log";
            xlabel('x')
            ylabel('KL(x)')
            legend([nse_h(1),nmem_h(1)],'Interpreter','latex')
            if save_figs
                saveas(bp, fullfile(fig_dir, [fig_name, '.png']))
            end

            fig_name = ['avg_kl_d_',convertStringsToChars(d_vec(d_idx)),'_s_',num2str(n_vec(n_idx))];
            figure('Name',fig_name, 'visible',fig_plot)
            hold on;
            nse_h = plot(quantiles, mean(kl_dist_nse_per_q,2), '-r', 'DisplayName', '$\hat{f}_i^{NSE}(x)$');
            nmem_h = plot(quantiles, mean(kl_dist_nmem_per_q,2), '-b', 'DisplayName', '$\hat{f}_i^{NMEM}(x)$');
            xlabel('x')
            ylabel('KL(x)')
            legend([nse_h(1),nmem_h(1)],'Interpreter','latex')
            bp = gca;
            if save_figs
                saveas(bp, fullfile(fig_dir, [fig_name, '.png']))
            end
        end

        % build table


        % MSE
        temp = vertcat(utils.reshape_groups(n_vec',mse_dist_nse),...
            utils.reshape_groups(n_vec',mse_dist_nmem));
    
        sample_power = temp(:,1);
        mse_total = temp(:,2);

        % KL
        temp = vertcat(utils.reshape_groups(n_vec',kl_dist_nse),...
            utils.reshape_groups(n_vec',kl_dist_nmem));
    
        kl_total = temp(:,2);
    
        % Distributions 
        distribution = repelem(d_vec(d_idx), length(temp(:,2)))';
        name = repelem(names(d_idx), length(temp(:,2)))';
    
        nse_label = repelem(["NSE"], size(utils.reshape_groups(n_vec',mse_dist_nse), 1));
        nmem_label = repelem(["NMEM"], size(utils.reshape_groups(n_vec',mse_dist_nmem), 1));


        estimator = vertcat(nse_label', nmem_label');

        
        dist_table = table(distribution, name, estimator, sample_power, mse_total,kl_total);
        dist_table = splitvars(dist_table);
    
        % append table per distribution to global table containing data 
        % for all distributions
        global_table = [global_table; dist_table];

    end

    % MSE per distribution
    fig_name = 'MSE_per_sample';
    figure('Name',fig_name, 'visible',fig_plot)
    b = boxchart(log(global_table.sample_power)/log(2), global_table.mse_total, 'GroupByColor',global_table.estimator);
    bp = gca;
    bp.XAxis.TickLabelInterpreter = 'latex';
    xlabel('$log_{2}(N)$','Interpreter','latex')
    ylabel('MSE','Interpreter','latex')
    legend('Location','northwest')
    if save_figs
        saveas(bp, fullfile(fig_dir, [fig_name, '.png']))
    end

    writetable(global_table,fullfile('data',table_name))

    for idx = 1:length(names)
        
        tester = global_table.name;
        mask = global_table.name == names(idx);
        data_dist = global_table(mask,:);
        
        fig_name = sprintf('mse_%s',names(idx));
        figure('Name',fig_name, 'visible',fig_plot)
        b = boxchart(log(data_dist.sample_power)/log(2), data_dist.mse_total, 'GroupByColor',data_dist.estimator);
        bp = gca;
        bp.XAxis.TickLabelInterpreter = 'latex';
        xlabel('$log_{2}(N)$','Interpreter','latex')
        ylabel('MSE','Interpreter','latex')
        legend('Location','northwest')
        title(convertCharsToStrings(names(idx)))
        if save_figs
            saveas(bp, fullfile(fig_dir, [fig_name, '.png']))
        end
    end

    % KL per distribution
    fig_name = 'KL_per_sample';
    figure('Name',fig_name, 'visible',fig_plot)
    b = boxchart(log(global_table.sample_power)/log(2), global_table.kl_total, 'GroupByColor',global_table.estimator);
    bp = gca;
    bp.XAxis.TickLabelInterpreter = 'latex';
    xlabel('$log_{2}(N)$','Interpreter','latex')
    ylabel('MSE','Interpreter','latex')
    legend('Location','northwest')
    if save_figs
        saveas(bp, fullfile(fig_dir, [fig_name, '.png']))
    end

    writetable(global_table,fullfile('data',table_name))

    for idx = 1:length(names)
        
        tester = global_table.name;
        mask = global_table.name == names(idx);
        data_dist = global_table(mask,:);
        
        fig_name = sprintf('kl_%s',names(idx));
        figure('Name',fig_name, 'visible',fig_plot)
        b = boxchart(log(data_dist.sample_power)/log(2), data_dist.kl_total, 'GroupByColor',data_dist.estimator);
        bp = gca;
        bp.XAxis.TickLabelInterpreter = 'latex';
        xlabel('$log_{2}(N)$','Interpreter','latex')
        ylabel('KL','Interpreter','latex')
        legend('Location','northwest')
        title(convertCharsToStrings(names(idx)))
        if save_figs
            saveas(bp, fullfile(fig_dir, [fig_name, '.png']))
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
            figure('Name',fig_name, 'visible',fig_plot)
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
                saveas(bp, fullfile(fig_dir, [fig_name, '.png']))
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
            figure('Name',fig_name, 'visible',fig_plot)
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
                saveas(bp, fullfile(fig_dir, [fig_name, '.png']))
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
            figure('Name',fig_name, 'visible',fig_plot)
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
                saveas(bp, fullfile(fig_dir, [fig_name, '.png']))
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
            figure('Name',fig_name, 'visible',fig_plot)
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
                saveas(bp, fullfile(fig_dir, [fig_name, '.png']))
            end
        end
    end
end




