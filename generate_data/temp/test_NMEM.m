close all;
clear all;
low = 1;
high = 1;
% 
% for i = low:high
% %     sample = importdata(sprintf(['sample_D_Uniform-Mix_T_','%d', '_S_8192.dat'],i));
% %     sample = importdata(sprintf(['sample_D_Uniform_T_','%d', '_S_1024.dat'],i));
%     sample = importdata(sprintf(['sample_D_Beta-a0p5-b0p5_T_','%d', '_S_1024.dat'],i));
%     
%     [failed, x, pdf, cdf, ~,lagrange] = EstimatePDF(sample);
%     solution{i}.x = x;
%     solution{i}.pdf = pdf;
% end

% sample = importdata(['D_Normal_T_3_S_4096_blockCDF_1.dat']);
sample = importdata(['sample_D_Beta-a0p5-b0p5_T_1_S_8192.dat']);

[sample, ind] = unique(sample);
% sample = sample(ind);

% [failed, x, pdf, cdf, u,sqr,lagrange] = EstimatePDF(sample);
[failed, x, pdf, cdf, sqr, lagrange, score, confidence, SURD] = EstimatePDF(sample);

% data = randn(1000, 1);
% x = min(data):0.1:max(data);
% f = normpdf(x);
% d = [x(:), f(:)];
% PDFAnalyze(data, 'distribution', d);

PDFAnalyze(sample, 'PlotType', 'sqr');

solution{1}.x = x;
solution{1}.pdf = pdf;
solution{1}.cdf = cdf;
solution{1}.sqr = sqr;

n = length(solution{1}.sqr);
dx = 1 / (n + 1);
u = dx:dx:(n * dx);
% idxOut = [find(sqr' > topThreshold), find(sqr' < bottomThreshold)];
% idxIn =  intersect(find(sqr' < topThreshold), find(sqr' > bottomThreshold));

solution{1}.u = u;

size(u)
size(sqr)

figure('Name','no transfrom NMEM pdfe')
hold on
for i = low:high
    plot(solution{i}.x,solution{i}.pdf, '-k')
    plot(solution{i}.x,solution{i}.cdf, '-r')
    plot(solution{i}.u,solution{i}.sqr, '-g')
end
xlim([0,1])
ylim([-1,1])
% ylim([0,6])

















