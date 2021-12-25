
low = 1;
high = 3;

for i = low:high
    sample = importdata(sprintf(['sample_D_Uniform-Mix_T_','%d', '_S_8192.dat'],i));
    sample = importdata(sprintf(['sample_D_Uniform_T_','%d', '_S_1024.dat'],i));
%     sample = importdata(sprintf(['sample_D_Beta-a0p5-b0p5_T_','%d', '_S_1024.dat'],i));
    
    [failed, x, pdf, cdf, ~,lagrange] = EstimatePDF(sample);
    solution{i}.x = x;
    solution{i}.pdf = pdf;
end

figure('Name','no transfrom NMEM pdfe')
hold on
for i = low:high
    plot(solution{i}.x,solution{i}.pdf, '-k')
end
xlim([0,10])
ylim([0,1])
% ylim([0,6])

















