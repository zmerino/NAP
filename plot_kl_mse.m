clc; clear; close all;

addpath("functions/")
addpath("functions_plotting/")

publicationQuality();

plot_mse_dist = false;
save_figs = true;

table_name = 'mse_kl_100.dat';
fig_dir = fullfile('figures_manuscript','kl_mse_quantiles');
write_dir = fullfile('data_3','kl_mse_data');

status = mkdir(write_dir);
status = mkdir(fig_dir);

data = readtable(fullfile(write_dir,table_name));
data_nap = data('NAP'==convertCharsToStrings(data.estimator),:);
data_nmem = data('NMEM'==convertCharsToStrings(data.estimator),:);

d_vec = ["Trimodal-Normal","Uniform","Normal","Beta-a0p5-b1p5","Beta-a2-b0p5","Beta-a0p5-b0p5","Generalized-Pareto","Stable"];
names = ["Trimodal-Normal","Uniform", "Normal", "Beta(0.5,1.5)", "Beta(2,0.5)", "Beta(0.5,0.5)", "Generalized-Pareto", "Stable"];
labels = {'$2^{10}$', '$2^{11}$', '$2^{12}$',...
    '$2^{13}$', '$2^{14}$', '$2^{15}$', '$2^{16}$', '$2^{17}$',...
    '$2^{18}$', '$2^{19}$', '$2^{20}$', '$2^{21}$', '$2^{22}$'};
n_vec = 2.^[8,9,10,11,12,13,14,15,16,17,18,19,20,21,22];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Quantiles %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% d_vec = ["Beta-a0p5-b1p5","Beta-a2-b0p5","Beta-a0p5-b0p5"];
trials = 100;
plot_q = false;

if plot_q
    for i = 1:length(d_vec)
        % save distribution name to object
        actual.dist_name = d_vec(i);

        for j = 1:length(n_vec)


            disp(['d: ', num2str(i), '/', num2str(length(d_vec)), ' n: ', num2str(j), '/', num2str(length(n_vec))])


            filename = sprintf('kl_per_q_%s_n_%s_t_%s.mat',d_vec(i), num2str(n_vec(j)), num2str(trials));

            load(fullfile(write_dir,['nse_', filename]), 'kl_dist_nse_per_q')
            load(fullfile(write_dir,['nmem_', filename]), 'kl_dist_nmem_per_q')

            filename = sprintf('mse_per_q_%s_n_%s_t_%s.mat',d_vec(i), num2str(n_vec(j)), num2str(trials));

            load(fullfile(write_dir,['nse_', filename]), 'mse_dist_nse_per_q')
            load(fullfile(write_dir,['nmem_', filename]), 'mse_dist_nmem_per_q')


            q_steps = size(kl_dist_nse_per_q,1);
            q_min = 0.01*100;
            q_max = 0.99*100;

            % KL ------------------------------

            nap_avg_kl = mean(kl_dist_nse_per_q,2);
            nmem_avg_kl = mean(kl_dist_nmem_per_q,2);
            qauntile = linspace(q_min,q_max,q_steps);

            fig_name = sprintf('kl_per_q_%s_n_%s_t_%s',d_vec(i), num2str(n_vec(j)), num2str(trials));
            figure('Name',fig_name)
            hold on;
            nse_h = plot(qauntile,nap_avg_kl, '-r', 'DisplayName', '$\hat{f}(Q^{-1}(y))_{NAP}$');
            nmem_h = plot(qauntile,nmem_avg_kl, '-b', 'DisplayName', '$\hat{f}(Q^{-1}(y))_{NMEM}$');
            legend([nse_h(1),nmem_h(1)],'Interpreter','latex')
            xlabel('$Quantile$', Interpreter='latex')
            ylabel('$ \langle KL \rangle $', Interpreter='latex')
            xtickformat('percentage')
            bp = gca;
            if save_figs
                saveas(bp, fullfile(fig_dir, [fig_name, '.png']))
            end


            % log scale
            fig_name = sprintf('log_kl_per_q_%s_n_%s_t_%s',d_vec(i), num2str(n_vec(j)), num2str(trials));
            figure('Name',fig_name)
            hold on;
            nse_h = plot(qauntile,nap_avg_kl, '-r', 'DisplayName', '$\hat{f}(Q^{-1}(y))_{NAP}$');
            nmem_h = plot(qauntile,nmem_avg_kl, '-b', 'DisplayName', '$\hat{f}(Q^{-1}(y))_{NMEM}$');
            legend([nse_h(1),nmem_h(1)],'Interpreter','latex')
            xlabel('$Quantile$', Interpreter='latex')
            ylabel('$ \langle KL \rangle $', Interpreter='latex')
            bp = gca;
            bp.YAxis.Scale ="log";
            if save_figs
                saveas(bp, fullfile(fig_dir, [fig_name, '.png']))
                saveas(bp, fullfile(fig_dir, [fig_name, '.fig']))
            end


            % MSE ------------------------------

            nap_avg_mse = mean(mse_dist_nse_per_q,2);
            nmem_avg_mse = mean(mse_dist_nmem_per_q,2);
            qauntile = linspace(q_min,q_max,q_steps);

            fig_name = sprintf('mse_per_q_%s_n_%s_t_%s',d_vec(i), num2str(n_vec(j)), num2str(trials));
            figure('Name',fig_name)
            hold on;
            nse_h = plot(qauntile,nap_avg_mse, '-r', 'DisplayName', '$\hat{f}(Q^{-1}(y))_{NAP}$');
            nmem_h = plot(qauntile,nmem_avg_mse, '-b', 'DisplayName', '$\hat{f}(Q^{-1}(y))_{NMEM}$');
            legend([nse_h(1),nmem_h(1)],'Interpreter','latex')
            xlabel('$Quantile$', Interpreter='latex')
            ylabel('$ \langle MSE \rangle $', Interpreter='latex')
            bp = gca;
            if save_figs
                saveas(bp, fullfile(fig_dir, [fig_name, '.png']))
                saveas(bp, fullfile(fig_dir, [fig_name, '.fig']))
            end


            % log scale
            fig_name = sprintf('log_mse_per_q_%s_n_%s_t_%s',d_vec(i), num2str(n_vec(j)), num2str(trials));
            figure('Name',fig_name)
            hold on;
            nse_h = plot(qauntile,nap_avg_mse, '-r', 'DisplayName', '$\hat{f}(Q^{-1}(y))_{NAP}$');
            nmem_h = plot(qauntile,nmem_avg_mse, '-b', 'DisplayName', '$\hat{f}(Q^{-1}(y))_{NMEM}$');
            legend([nse_h(1),nmem_h(1)],'Interpreter','latex')
            xlabel('$Quantile$', Interpreter='latex')
            ylabel('$ \langle MSE \rangle $', Interpreter='latex')
            bp = gca;
            bp.YAxis.Scale ="log";
            if save_figs
                saveas(bp, fullfile(fig_dir, [fig_name, '.png']))
                saveas(bp, fullfile(fig_dir, [fig_name, '.fig']))
            end

        end

    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MSE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fig_name = 'MSE All Distributions';
figure('Name',fig_name)
b = boxchart(log(data.sample_power)/log(2), data.mse_total, 'GroupByColor',data.estimator);
bp = gca;
bp.XAxis.TickLabelInterpreter = 'latex';
xlabel('$log_{2}(N)$','Interpreter','latex')
ylabel('MSE','Interpreter','latex')
legend('Location','northwest')
if save_figs
    saveas(bp, fullfile(fig_dir, [fig_name, '.png']))
    saveas(bp, fullfile(fig_dir, [fig_name, '.fig']))
end

fig_name = 'NAP MSE Per Distribution';
figure('Name',fig_name)
b = boxchart(log(data_nap.sample_power)/log(2), data_nap.mse_total, 'GroupByColor',data_nap.distribution);
bp = gca;
bp.XAxis.TickLabelInterpreter = 'latex';
xlabel('$log_{2}(N)$','Interpreter','latex')
ylabel('MSE','Interpreter','latex')
legend('Location','northwest')
if save_figs
    saveas(bp, fullfile(fig_dir, [fig_name, '.png']))
    saveas(bp, fullfile(fig_dir, [fig_name, '.fig']))
end

fig_name = 'NMEM MSE Per Distribution';
figure('Name',fig_name)
b = boxchart(log(data_nmem.sample_power)/log(2), data_nmem.mse_total, 'GroupByColor',data_nmem.distribution);
bp = gca;
bp.XAxis.TickLabelInterpreter = 'latex';
xlabel('$log_{2}(N)$','Interpreter','latex')
ylabel('MSE','Interpreter','latex')
legend('Location','northwest')
if save_figs
    saveas(bp, fullfile(fig_dir, [fig_name, '.png']))
    saveas(bp, fullfile(fig_dir, [fig_name, '.fig']))
end



for idx = 1:length(names)

    tester = data.name;
    mask = data.name == names(idx);
    data_dist = data(mask,:);
    data_dist.kl_total = log10(data_dist.kl_total);

    fig_name = sprintf('mse_%s',names(idx));
    figure('Name',fig_name)
    b = boxchart(log(data_dist.sample_power)/log(2), data_dist.mse_total, 'GroupByColor',data_dist.estimator);
    bp = gca;
    bp.XAxis.TickLabelInterpreter = 'latex';
    xlabel('$log_{2}(N)$','Interpreter','latex')
    ylabel('KL','Interpreter','latex')
    legend('Location','northwest')
    title(convertCharsToStrings(names(idx)))
    %     bp.YAxis.Scale ="log";
    if save_figs
        saveas(bp, fullfile(fig_dir, [fig_name, '.png']))
        saveas(bp, fullfile(fig_dir, [fig_name, '.fig']))
    end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% KL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fig_name = 'KL All Distributions';
figure('Name',fig_name)
b = boxchart(log(data.sample_power)/log(2), data.kl_total, 'GroupByColor',data.estimator);
bp = gca;
bp.XAxis.TickLabelInterpreter = 'latex';
xlabel('$log_{2}(N)$','Interpreter','latex')
ylabel('KL','Interpreter','latex')
legend('Location','northwest')
if save_figs
    saveas(bp, fullfile(fig_dir, [fig_name, '.png']))
    saveas(bp, fullfile(fig_dir, [fig_name, '.fig']))
end


fig_name = 'NAP KL Per Distribution';
figure('Name',fig_name)
b = boxchart(log(data_nap.sample_power)/log(2), data_nap.kl_total, 'GroupByColor',data_nap.distribution);
bp = gca;
bp.XAxis.TickLabelInterpreter = 'latex';
xlabel('$log_{2}(N)$','Interpreter','latex')
ylabel('KL','Interpreter','latex')
legend('Location','northwest')
if save_figs
    saveas(bp, fullfile(fig_dir, [fig_name, '.png']))
    saveas(bp, fullfile(fig_dir, [fig_name, '.fig']))
end

fig_name = 'NMEM KL Per Distribution';
figure('Name',fig_name)
b = boxchart(log(data_nmem.sample_power)/log(2), data_nmem.kl_total, 'GroupByColor',data_nmem.distribution);
bp = gca;
bp.XAxis.TickLabelInterpreter = 'latex';
xlabel('$log_{2}(N)$','Interpreter','latex')
ylabel('KL','Interpreter','latex')
legend('Location','northwest')
if save_figs
    saveas(bp, fullfile(fig_dir, [fig_name, '.png']))
    saveas(bp, fullfile(fig_dir, [fig_name, '.fig']))
end

for idx = 1:length(names)

    tester = data.name;
    mask = data.name == names(idx);
    data_dist = data(mask,:);
%     data_dist.kl_total = log10(data_dist.kl_total);

    fig_name = sprintf('kl_%s',names(idx));
    figure('Name',fig_name)
    b = boxchart(log(data_dist.sample_power)/log(2), data_dist.kl_total, 'GroupByColor',data_dist.estimator);
    bp = gca;
    bp.XAxis.TickLabelInterpreter = 'latex';
    xlabel('$log_{2}(N)$','Interpreter','latex')
    ylabel('KL','Interpreter','latex')
    legend('Location','northwest')
    title(convertCharsToStrings(names(idx)))
    %     bp.YAxis.Scale ="log";
    if save_figs
        saveas(bp, fullfile(fig_dir, [fig_name, '.png']))
        saveas(bp, fullfile(fig_dir, [fig_name, '.fig']))
    end
end






