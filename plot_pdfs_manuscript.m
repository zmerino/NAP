clc; clear; close all;

addpath("functions/")
addpath("functions_plotting/")



%%%%%%%%%%%%%%%%%%%%%%%%%%%%% To Plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

plot_mse_dist = true;
table_name = 'mse_kl_set0.dat';

plot_pdf = false;
plot_cdf = false;
plot_sqr = false;
plot_heavy = false;
save_figs = true;

% figure directory
fig_dir = fullfile('figures_manuscript','100_trials_v2');
status = mkdir(fig_dir);

% choose to visualise figures or not
fig_plot = 'on';

%%%%%%%%%%%%%%%%%%%%%%%%%%% Import data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Create actual distrobution onjects
actual = distributions;

% Path to directory
% dir_name = fullfile('data','pdf_estimates');
dir_name = fullfile('data','cpu_20_t_50_maxN_22_pdf_estimates');

dir_name = fullfile('data','cpu_20_t_50_set_1');
dir_name = fullfile('data','cpu_20_t_100');

% Define the etimates to plot

% Sample range
% n_vec = 2.^[8,9,10,11,12,13,14,15,16,17,18];
n_vec = 2.^[14];
n_vec = 2.^[8,9,10,11,12,13,14,15,16,17,18,19,20,21,22];
n_vec = 2.^[10,15];
n_vec = 2.^[20];
n_vec = 2.^[18];
% Distribution range
% d_vec = ["Beta-a0p5-b1p5","Beta-a2-b0p5","Beta-a0p5-b0p5","Generalized-Pareto","Stable"];
% d_vec = ["Trimodal-Normal","Uniform","Normal","Uniform-Mix","Beta-a0p5-b1p5","Beta-a2-b0p5","Beta-a0p5-b0p5","Generalized-Pareto","Stable"];
d_vec = ["Trimodal-Normal","Uniform","Normal","Beta-a0p5-b1p5","Beta-a2-b0p5","Beta-a0p5-b0p5","Generalized-Pareto","Stable"];
% d_vec = ["Trimodal-Normal","Uniform","Normal"];
% d_vec = ["Beta-a0p5-b1p5","Beta-a2-b0p5","Beta-a0p5-b0p5","Generalized-Pareto"];
% d_vec = ["Normal"];

% Define labels for figures
% names = {'Uniform-Mix', 'Generalized-Pareto', 'Stable',...
%     'Tri-Modal-Normal', 'Normal', 'Uniform', 'Beta(0.5,1.5)',...
%     'Beta(2,0.5)', 'Beta(0.5,0.5)'}';
% names = {'Beta(0.5,1.5)','Beta(2,0.5)', 'Beta(0.5,0.5)', 'Generalized-Pareto', 'Stable'}';
names = ["Trimodal-Normal","Uniform","Normal","Beta(0.5,1.5)","Beta(2,0.5)", "Beta(0.5,0.5)", "Generalized-Pareto","Stable"];

d_vec = ["Trimodal-Normal","Uniform","Beta-a0p5-b1p5","Beta-a2-b0p5","Beta-a0p5-b0p5","Generalized-Pareto","Stable"];
names = ["Trimodal-Normal","Uniform","Beta(0.5,1.5)","Beta(2,0.5)", "Beta(0.5,0.5)", "Generalized-Pareto","Stable"];
% d_vec = ["Beta-a0p5-b0p5","Generalized-Pareto","Stable"];
% names = ["Beta(0.5,0.5)", "Generalized-Pareto","Stable"];
% 
% % heavy tails
% d_vec = ["Generalized-Pareto","Stable"];
% names = ["Generalized-Pareto","Stable"];
% d_vec = ["Trimodal-Normal","Normal"];
% names = ["Trimodal-Normal","Normal"];


% 
% d_vec = ["Beta-a0p5-b1p5","Beta-a2-b0p5","Beta-a0p5-b0p5"];
% names = ["Uniform","Normal","Beta(0.5,1.5)","Beta(2,0.5)", "Beta(0.5,0.5)"];

% Trials per sample
trials = 100;
plt_trials = 10;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Load data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialize empty structure to store data in
n = length(n_vec);
d = length(d_vec);
t = length(d_vec);
nap = cell(d,n,t,7);
nmem = cell(d,n,t,7);

for i = 1:length(d_vec)
    for j = 1:length(n_vec)
        for k = 1:trials
                        
            % Make sure that distribution names are saved as character vectors

            % track calculations
            disp(['dist: ',num2str(i),'/',num2str(length(d_vec)),...
                ' s: ',num2str(j),'/',num2str(length(n_vec)),...
                ' t: ',num2str(k),'/',num2str(trials)])

            % distribution names
            nap{i,j,k,1} = d_vec(i);
            nmem{i,j,k,1} = d_vec(i);
            % sample size
            nap{i,j,k,2} = n_vec(j);
            nmem{i,j,k,2} = n_vec(j);

            % Load data from estimate directory
            load(fullfile(dir_name,['nse_',num2str(d_vec(i)),'_t_',num2str(trials),'_s_', num2str(n_vec(j)),'.mat']));
            load(fullfile(dir_name,['nmem_',num2str(d_vec(i)),'_t_',num2str(trials),'_s_', num2str(n_vec(j)),'.mat']));

            % x
            nse_data = nse_pdf_data;
            nmem_data = nmem_pdf_data;

            nap{i,j,k,3} = nse_data{2,k}'; % make sure column vector
            nmem{i,j,k,3} = nmem_data{2,k};

            % pdf
            nap{i,j,k,4} = nse_data{3,k}'; % make sure column vector
            nmem{i,j,k,4} = nmem_data{3,k};

            % cdf
            nap{i,j,k,5} = utils_analysis.trapz(nse_data{2,k}, nse_data{3,k});
            nmem{i,j,k,5} = utils_analysis.trapz(nmem_data{2,k}, nmem_data{3,k});

            % u and sqr
            [nap{i,j,k,6}, nap{i,j,k,7}] = utils.sqr(nse_data{2,k},nse_data{3,k},nse_data{1,k});
            [nmem{i,j,k,6}, nmem{i,j,k,7}] = utils.sqr(nmem_data{2,k},nmem_data{3,k},nmem_data{1,k});

        end
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%% Ploting data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% MSE Plots %%%
actual = distributions;
actual.generate_data = false;
% tracks kl/utils_analysis.mse for entire sample
global_table = table();

if plot_mse_dist
    
    qnum = 50;
%     qnum = 20;
%     qnum = 100;

    % x values for given quantiles
    for d_idx = 1:size(nap, 1)

        % define distribution name for plotting and building a data table
        actual.dist_name = d_vec(d_idx);

        mse_dist_nmem = zeros(size(nap, 2),plt_trials);
        kl_dist_nmem = zeros(size(nap, 2),plt_trials);
        mse_dist_nse = zeros(size(nap, 2),plt_trials);
        kl_dist_nse = zeros(size(nap, 2),plt_trials);
        for n_idx = 1:size(nap, 2)

            mse_dist_nse_per_q = zeros(qnum,plt_trials);
            mse_dist_nmem_per_q = zeros(qnum,plt_trials);
            kl_dist_nmem_per_q = zeros(qnum,plt_trials);
            kl_dist_nse_per_q = zeros(qnum,plt_trials);
            for t_idx = 1:plt_trials

                % track calculations
                disp(['quantile dist: ',num2str(d_idx),'/',num2str(size(nap, 1)),...
                    ' s: ',num2str(n_idx),'/',num2str(size(nap, 2)),...
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


                if sum(~isfinite(xs)) || sum(~isfinite(fs))|| sum(~isfinite(Fs))

                    warning('non-finite values')

                    test = 'test';
                end
                
               % MSE per quantile
                [mse_dist_nmem_per_q(:,t_idx), quantiles,xq,fq] = utils_analysis.quant_metric(xs, fs, Fs, qnum, d_vec(d_idx),'MSE');
                % MSE per trial
                mse_dist_nmem(n_idx, t_idx) = utils_analysis.mse(actual.pdf_y, fs);
                % KL per quantile
                [kl_dist_nmem_per_q(:,t_idx), ~,~,~] = utils_analysis.quant_metric(xs, fs, Fs, qnum, d_vec(d_idx),'KL');
                % KL per trial
                kl_dist_nmem(n_idx, t_idx) = utils_analysis.kl(actual.pdf_y, fs');

               % NAP ------------------------------------------------------

                xs = nap{d_idx,n_idx,t_idx,3}(:,1);
                fs = nap{d_idx,n_idx,t_idx,4}(:,1);
                Fs = nap{d_idx,n_idx,t_idx,5}(:,1);
                % get actual distribution
                actual.min_limit = min(xs);
                actual.max_limit = max(xs);
                actual.x = xs;
                actual = actual.dist_list();

               % MSE per quantile
                [mse_dist_nse_per_q(:,t_idx),~,~,~] = utils_analysis.quant_metric(xs, fs, Fs, qnum, d_vec(d_idx),'MSE');
                % MSE per trial
                mse_dist_nse(n_idx, t_idx) = utils_analysis.mse(actual.pdf_y, fs);
                % KL per quantile
                [kl_dist_nse_per_q(:,t_idx), ~,~,~] = utils_analysis.quant_metric(xs, fs, Fs, qnum, d_vec(d_idx),'KL');
                % KL per trial
                kl_dist_nse(n_idx, t_idx) = utils_analysis.kl(actual.pdf_y, fs');


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
                saveas(bp, fullfile(fig_dir, [fig_name, '.fig']))
            end

            fig_name = ['mse_d_',convertStringsToChars(d_vec(d_idx)),'_s_',num2str(n_vec(n_idx))];
            figure('Name',fig_name, 'visible',fig_plot)
            hold on;
            for t = 1:size(mse_dist_nse_per_q,2)
                nse_h = plot(quantiles, mse_dist_nse_per_q(:,t), '-r', 'DisplayName', '$\hat{f}_i^{NAP}(x)$');
            end
            for t = 1:size(mse_dist_nse_per_q,2)
                nmem_h = plot(quantiles, mse_dist_nmem_per_q(:,t), '-b', 'DisplayName', '$\hat{f}_i^{NMEM}(x)$');
            end
            bp = gca;
%             bp.YAxis.Scale ="log";
            xlabel('x')
            ylabel('MSE(x)')
            legend([nse_h(1),nmem_h(1)],'Interpreter','latex')
            if save_figs
                saveas(bp, fullfile(fig_dir, [fig_name, '.png']))
                saveas(bp, fullfile(fig_dir, [fig_name, '.fig']))
            end


            fig_name = ['avg_mse_d_',convertStringsToChars(d_vec(d_idx)),'_s_',num2str(n_vec(n_idx))];
            figure('Name',fig_name, 'visible',fig_plot)
            hold on;
            nse_h = plot(quantiles, mean(mse_dist_nse_per_q,2), '-r', 'DisplayName', '$\hat{f}_i^{NAP}(x)$');
            nmem_h = plot(quantiles, mean(mse_dist_nmem_per_q,2), '-b', 'DisplayName', '$\hat{f}_i^{NMEM}(x)$');
            xlabel('x')
            ylabel('MSE(x)')
            legend([nse_h(1),nmem_h(1)],'Interpreter','latex')
            bp = gca;
%             bp.YAxis.Scale ="log";
            if save_figs
                saveas(bp, fullfile(fig_dir, [fig_name, '.png']))
                saveas(bp, fullfile(fig_dir, [fig_name, '.fig']))
            end

            fig_name = ['kl_d_',convertStringsToChars(d_vec(d_idx)),'_s_',num2str(n_vec(n_idx))];
            figure('Name',fig_name, 'visible',fig_plot)
            hold on;
            for t = 1:size(kl_dist_nmem_per_q,2)
                nse_h = plot(quantiles, kl_dist_nse_per_q(:,t), '-r', 'DisplayName', '$\hat{f}_i^{NAP}(x)$');
                nmem_h = plot(quantiles, kl_dist_nmem_per_q(:,t), '-b', 'DisplayName', '$\hat{f}_i^{NMEM}(x)$');
            end
            bp = gca;
%             bp.YAxis.Scale ="log";
            xlabel('x')
            ylabel('KL(x)')
            legend([nse_h(1),nmem_h(1)],'Interpreter','latex')
            if save_figs
                saveas(bp, fullfile(fig_dir, [fig_name, '.png']))
                saveas(bp, fullfile(fig_dir, [fig_name, '.fig']))
            end

            fig_name = ['avg_kl_d_',convertStringsToChars(d_vec(d_idx)),'_s_',num2str(n_vec(n_idx))];
            figure('Name',fig_name, 'visible',fig_plot)
            hold on;
            nse_h = plot(quantiles, mean(kl_dist_nse_per_q,2), '-r', 'DisplayName', '$\hat{f}_i^{NAP}(x)$');
            nmem_h = plot(quantiles, mean(kl_dist_nmem_per_q,2), '-b', 'DisplayName', '$\hat{f}_i^{NMEM}(x)$');
            xlabel('x')
            ylabel('KL(x)')
            legend([nse_h(1),nmem_h(1)],'Interpreter','latex')
            bp = gca;
%             bp.YAxis.Scale ="log";
            if save_figs
                saveas(bp, fullfile(fig_dir, [fig_name, '.png']))
                saveas(bp, fullfile(fig_dir, [fig_name, '.fig']))
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
    
        nse_label = repelem(["NAP"], size(utils.reshape_groups(n_vec',mse_dist_nse), 1));
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
        ylabel('CPU Time [$log_{10}(sec)$]','Interpreter','latex')
        legend('Location','northwest')
        title(convertCharsToStrings(names(idx)))
        if save_figs
            saveas(bp, fullfile(dir_name, [fig_name, '.png']))
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
        ylabel('CPU Time [$log_{10}(sec)$]','Interpreter','latex')
        legend('Location','northwest')
        title(convertCharsToStrings(names(idx)))
        if save_figs
            saveas(bp, fullfile(dir_name, [fig_name, '.png']))
        end
    end

end


%%% PDFs %%%
if plot_pdf
    % Plot data the order of indices corresponds to distribution, sample size,
    % trial. Plot all trials per sample per distrobution
    for d_idx = 1:size(nap, 1)

        for n_idx = 1:size(nap, 2)

            fig_name = ['nmem_PDF_d',convertStringsToChars(d_vec(d_idx)),num2str(d_idx),'_s_',num2str(n_vec(n_idx)),'_t_',num2str(trials)];
            figure('Name',fig_name, 'visible',fig_plot)
            hold on;
            % Set actual distribution type
            actual.dist_name = nap{d_idx,n_idx,1,1}(1,:);
            % Set distrobution over specific range i.e. the x values from the
            % first trial estimate

            actual.x = nap{d_idx,n_idx,1,3}(:,1);
            % Generate distrobution attributes
            actual = dist_list(actual);

            % Remove zeros in actual.pdf_y generated from inf/nan values
            [row,~] = find(actual.pdf_y > 0);
            for s_idx = 1:trials
                nmem_h = plot(nmem{d_idx,n_idx,s_idx,3}, nmem{d_idx,n_idx,s_idx,4},'-', 'DisplayName', '$\hat{f}_i^{NMEM}(x)$');
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
            legend([g],'Interpreter','latex')
            if save_figs
                saveas(bp, fullfile(fig_dir, [fig_name, '.png']))
            end


            fig_name = ['nap_PDF_d',convertStringsToChars(d_vec(d_idx)),num2str(d_idx),'_s_',num2str(n_vec(n_idx)),'_t_',num2str(trials)];
            figure('Name',fig_name, 'visible',fig_plot)
            hold on;
            % Set actual distribution type
            actual.dist_name = nap{d_idx,n_idx,1,1}(1,:);
            % Set distrobution over specific range i.e. the x values from the
            % first trial estimate

            actual.x = nap{d_idx,n_idx,1,3}(:,1);
            % Generate distrobution attributes
            actual = dist_list(actual);

            % Remove zeros in actual.pdf_y generated from inf/nan values
            [row,~] = find(actual.pdf_y > 0);
            for s_idx = 1:trials
                nse_h = plot(nap{d_idx,n_idx,s_idx,3}, nap{d_idx,n_idx,s_idx,4},'-', 'DisplayName', '$\hat{f}_i^{NAP}(x)$');
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
            legend([g],'Interpreter','latex')
            if save_figs
                saveas(bp, fullfile(fig_dir, [fig_name, '.png']))
            end


            fig_name = ['combined_PDF_d',convertStringsToChars(d_vec(d_idx)),num2str(d_idx),'_s_',num2str(n_vec(n_idx)),'_t_',num2str(trials)];
            figure('Name',fig_name, 'visible',fig_plot)
            hold on;
            % Set actual distribution type
            actual.dist_name = nap{d_idx,n_idx,1,1}(1,:);
            % Set distrobution over specific range i.e. the x values from the
            % first trial estimate

            actual.x = nap{d_idx,n_idx,1,3}(:,1);
            % Generate distrobution attributes
            actual = dist_list(actual);

            % Remove zeros in actual.pdf_y generated from inf/nan values
            [row,~] = find(actual.pdf_y > 0);
            for s_idx = 1:trials
                nmem_h = plot(nmem{d_idx,n_idx,s_idx,3}, nmem{d_idx,n_idx,s_idx,4},'-b', 'DisplayName', '$\hat{f}_i^{NMEM}(x)$');
            end
            for s_idx = 1:trials
                nse_h = plot(nap{d_idx,n_idx,s_idx,3}, nap{d_idx,n_idx,s_idx,4},'-r', 'DisplayName', '$\hat{f}_i^{NAP}(x)$');
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
    for d_idx = 1:size(nap, 1)

        for n_idx = 1:size(nap, 2)

            fig_name = ['CDF_d',convertStringsToChars(d_vec(d_idx)),num2str(d_idx),'_s_',num2str(n_vec(n_idx)),'_t_',num2str(trials)];
            figure('Name',fig_name, 'visible',fig_plot)
            hold on;
            % Set actual distribution type
            actual.dist_name = nap{d_idx,n_idx,1,1};
            % Set distrobution over specific range i.e. the x values from the
            % first trial estimate
            actual.x = nap{d_idx,n_idx,1,3}(:,1);
            % Generate distrobution attributes
            actual = dist_list(actual);

            % Remove zeros in actual.pdf_y generated from inf/nan values
            [row,~] = find(actual.pdf_y > 0);

            for s_idx = 1:trials
                nmem_h = plot(nmem{d_idx,n_idx,s_idx,3}, nmem{d_idx,n_idx,s_idx,5},'-b', 'DisplayName', '$\hat{F}_i^{NMEM}(x)$');
            end
            for s_idx = 1:trials
                nse_h = plot(nap{d_idx,n_idx,s_idx,3}, nap{d_idx,n_idx,s_idx,5},'-r', 'DisplayName', '$\hat{F}_i^{NAP}(x)$');
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
    for d_idx = 1:size(nap, 1)

        for n_idx = 1:size(nap, 2)

            fig_name = ['nmem_SQR_d',convertStringsToChars(d_vec(d_idx)),num2str(d_idx),'_s_',num2str(n_vec(n_idx)),'_t_',num2str(trials)];
            figure('Name',fig_name, 'visible',fig_plot)
            hold on;
            lemondrop();
            % Set actual distribution type
            actual.dist_name = nap{d_idx,n_idx,1,1};
            % Set distrobution over specific range i.e. the x values from the
            % first trial estimate
            actual.x = nap{d_idx,n_idx,1,3}(:,1);
            % Generate distrobution attributes
            actual = dist_list(actual);

            % Remove zeros in actual.pdf_y generated from inf/nan values
            [row,~] = find(actual.pdf_y > 0);

            for s_idx = 1:trials
                nmem_h = plot(nmem{d_idx,n_idx,s_idx,6}, nmem{d_idx,n_idx,s_idx,7},'-', 'DisplayName', '$\hat{SQR}_i^{NMEM}(x)$');
            end
            bp = gca;
            xlabel('$x$','Interpreter','latex')
            ylabel('$\hat{SQR}(x)$','Interpreter','latex')
%             legend([nmem_h(1)],'Interpreter','latex')
            xlabel('$x$','Interpreter','latex')
            ylabel('$sqr(x)$','Interpreter','latex')
            if save_figs
                saveas(bp, fullfile(fig_dir, [fig_name, '.png']))
            end


            fig_name = ['nap_SQR_d',convertStringsToChars(d_vec(d_idx)),num2str(d_idx),'_s_',num2str(n_vec(n_idx)),'_t_',num2str(trials)];
            figure('Name',fig_name, 'visible',fig_plot)
            hold on;
            lemondrop();
            % Set actual distribution type
            actual.dist_name = nap{d_idx,n_idx,1,1};
            % Set distrobution over specific range i.e. the x values from the
            % first trial estimate
            actual.x = nap{d_idx,n_idx,1,3}(:,1);
            % Generate distrobution attributes
            actual = dist_list(actual);

            % Remove zeros in actual.pdf_y generated from inf/nan values
            [row,~] = find(actual.pdf_y > 0);

            for s_idx = 1:trials
                nse_h = plot(nap{d_idx,n_idx,s_idx,6}, nap{d_idx,n_idx,s_idx,7},'-', 'DisplayName', '$\hat{SQR}_i^{NAP}(x)$');
            end
            bp = gca;
            xlabel('$x$','Interpreter','latex')
            ylabel('$\hat{SQR}(x)$','Interpreter','latex')
%             legend([nse_h(1)],'Interpreter','latex')
            xlabel('$x$','Interpreter','latex')
            ylabel('$SQR(x)$','Interpreter','latex')
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
    for d_idx = 1:size(nap, 1)

        for n_idx = 1:size(nap, 2)

            fig_name = ['Heavy_Tails_d',convertStringsToChars(d_vec(d_idx)),num2str(d_idx),'_s_',num2str(n_vec(n_idx)),'_t_',num2str(trials)];
            figure('Name',fig_name, 'visible',fig_plot)
            hold on;
            for j = 1:size(nap{d_idx,n_idx,3}, 2)

                % Set actual distribution type
                actual.dist_name = nap{d_idx,n_idx,1,1};
%                 % Set distrobution over specific range i.e. the x values from the
%                 % first trial estimate
%                 xs = nap{d_idx,n_idx,1,3}(:,1);
%                 x_min = min(xs);
%                 x_max = max(xs);
% 
%                 actual.x = linspace(x_min,x_max*2,1000);
%                 % Generate distrobution attributes
%                 actual = dist_list(actual);

                % Remove zeros in actual.pdf_y generated from inf/nan values
%                 [row,~] = find(actual.pdf_y > 0);

                for s_idx = 1:trials
                    % Set distrobution over specific range i.e. the x values from the
                    % first trial estimate
                    xs = nap{d_idx,n_idx,1,3}(:,1);
                    x_min = min(xs);
                    x_max = max(xs);

                    actual.x = linspace(x_min,x_max*2,1000);
                    % Generate distrobution attributes
                    actual = dist_list(actual);

                    [nse_htx, nse_hty, nse_htx_act, nse_hty_act] = utils_analysis.heavy(nap{d_idx,n_idx,s_idx,3}(:,j),...
                        nap{d_idx,n_idx,s_idx,5}(:,j), actual);

                    [nmem_htx, nmem_hty, nmem_htx_act, nmem_hty_act] = utils_analysis.heavy(nmem{d_idx,n_idx,s_idx,3}(:,j),...
                        nmem{d_idx,n_idx,s_idx,5}(:,j), actual);

                    nse_h = plot(real(nse_htx), real(nse_hty), '-r', 'DisplayName', '$\hat{F}_i^{NAP}(x)$');
                    nmem_h = plot(real(nmem_htx), real(nmem_hty), '-b', 'DisplayName', '$\hat{F}_i^{NMEM}(x)$');

                    xs = nap{d_idx,n_idx,s_idx,3}(:,1);

                    % update x-limits for actual cdf
                    if s_idx == 1
                        x_min_act = min(nse_htx_act);
                        x_max_act = max(nse_htx_act);
                    else
                        if min(nse_htx_act) < x_min_act
                            x_min_act = min(nse_htx_act);
                        end
                        if max(nse_htx_act) > x_max_act
                            x_max_act = max(nse_htx_act);
                        end
                    end

                end

                actual.x = linspace(x_min_act*2,x_max_act*2,1000);
                % Generate distrobution attributes
                actual = dist_list(actual);


                htx_act = log(actual.x);
                hty_act = log(1-actual.cdf_y);

                g = plot(nse_htx_act,nse_hty_act,'-k', 'DisplayName','$F(x)$');
            end
            bp = gca;
            set(bp,  'YScale', 'log')
            xlabel('$x$','Interpreter','latex')
            ylabel('$log(1 - \hat{F}(x))$','Interpreter','latex')
            legend([nse_h(1),nmem_h(1), g],'Interpreter','latex', 'Location','southwest')
            if save_figs
                saveas(bp, fullfile(fig_dir, [fig_name, '.png']))
            end
        end
    end
end




