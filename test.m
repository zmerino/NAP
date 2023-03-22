



% distribution_vector = ["Normal","Beta-a0p5-b0p5"];
% names = ["Normal", "Beta(0.5,0.5)"];

distribution_vector = ["Beta-a0p5-b0p5"];
names = ["Beta(0.5,0.5)"];
trials = '3';
min_pow = '13';
max_pow = '17';
cpu_n = '3';
cpp_code = "cpp_code/";


% 
% driver_func(distribution_vector,names,trials, min_pow,max_pow,cpu_n,cpp_code)
driver_func_n(distribution_vector,names,trials, min_pow,max_pow,cpu_n,cpp_code)


% distribution_vector = ["Normal","Beta-a0p5-b0p5"];
% names = ["Normal", "Beta(0.5,0.5)"];
% 
% trials = '2';
% min_pow = '13';
% max_pow = '14';
% cpu_n = '20';
% cpp_code = "cpp_code/";
% 
% for i =1:length(distribution_vector)
%     driver_func(distribution_vector(i),names(i),trials, min_pow,max_pow,cpu_n,cpp_code)
% end