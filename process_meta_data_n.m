
close all;

addpath("functions/")

save_figs = true;

% Import data table
file_dir = fullfile('data_large_n');
write_data = fullfile('data_large_n', 'meta_data_1.dat');

% distribution_vector = ["Trimodal-Normal","Uniform","Normal","Beta-a0p5-b1p5","Beta-a2-b0p5","Beta-a0p5-b0p5","Generalized-Pareto","Stable","Stable2","Stable3"];
distribution_vector = ["Beta-a0p5-b1p5","Generalized-Pareto"];
% distribution_vector = ["Beta-a0p5-b1p5","Beta-a2-b0p5","Beta-a0p5-b0p5"];
trials = 1;
cpu_n = 40;

large_table = table();

for idx = 1:length(distribution_vector)

    disp('PROCESSING DATA TABLE:')
    % construct data table name per distribution
    filename = sprintf('%s_cpu_%s_t_%s.dat',distribution_vector(idx), num2str(cpu_n), num2str(trials))
    % file path to data table
    file_path = fullfile(file_dir,filename);

    % load data table for a given distribution
    data = readtable(file_path);

    % append data table to global data table
    large_table = [large_table; data];

end


writetable(large_table, write_data)