% Master driver to test different subsample sizes and total number of
% subsamples


percSample = [0.2,0.4,0.6,0.8,0.9];
numSubs = [10,20,40,60,80,100];


CT = zeros(length(numSubs),length(percSample));

for i = 1:length(percSample)
     
    cpuHold = [];
    for j = 1:length(numSubs)
        
        tintial = cputime;
        stitchPDF_driver(percSample(i),numSubs(j))
        tcpu = cputime-tintial;
        
        cpuHold = [cpuHold,tcpu];
        
    end
    
    CT(:,i) = cpuHold;
       
end

dlmwrite('Subsample_times.txt',CT,'delimiter',' ');

figure('Name','Optimum parameters')
surf(percSample,numSubs,CT)
pngfile = strcat('optimum_plot.png');
saveas(gcf,pngfile)
figfile = strcat('optimum_plot.fig');
saveas(gcf,figfile)