% function stitchPDF_Driver(percSample,numSubs)
% ^^ for use with Master_driver.m to different subsampling parameters

% track command window
diary commandWindowOUT.txt

% use when stitchPDF_Driver.m is a script and not a function
clc;clear all; close all;

% % % % % % % % % % % % %
% % % percSample = 1;   %
% % % numSubs = 1;      %
% % % % % % % % % % % % %
% use ^^^ when one wishes to bypass subsampling routine

% SUBSAMPLING PARAMETERS---------------------------------------------------
% percentage of sample used to create subsample
percSample = 1;
% number of subsamples to generate
numSubs = 1;
%--------------------------------------------------------------------------

tic
%% User Options

% Script switching board %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SEcallflag =                true;   %<- true/false call SE on/off
SEplotflag =                true;   %<- true/false plot SE results on/off
dataTypeflag =              true;   %<- true/false integer powers of 2/real powers of 2
savePNG =                   true;   %<- true/false save .png of plots on/off
% random data generation parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
maxSamplesExp =             20; %<---- maximum exponent to generate samples
minSamplesExp =             10; %<---- minimum exponent to generate samples
precision =                 15; %<---- condtrol number of digits for created data
% will need to update code to handle trials and keep SRMSE
ntrials =                   1;  %<--- trials to run to generate heuristics for programs
step =                      1;  %<---- control synthetic random samples to skip being created
lowLim =                    0;  %<--- lower limit to plot
upLim =                     10; %<--- upper limit to plot
% Example distribution to test %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % distributionVector = ["Square-periodic","Uniform-Mix","Normal-Contaminated","Beta-a0p5-b0p5","Beta-a0p5-b1p5","Beta-a2-b0p5",...
% %     "Bimodal-Normal","Trimodal-Normal","BirnbaumSaunders"];
% distributionVector = ["Uniform-Mix","Normal-Contaminated","Bimodal-Normal","Trimodal-Normal",...
%     "BirnbaumSaunders","Beta-a0p5-b0p5","Beta-a0p5-b1p5","Beta-a2-b0p5"];
% distributionVector = ["Beta-a0p5-b1p5","Beta-a2-b0p5","Beta-a0p5-b0p5",...
%     "Bimodal-Normal","BirnbaumSaunders","Burr",...
%     "Exponential","Extreme-Value","Gamma","Generalized-Extreme-Value",...
%     "Generalized-Pareto","HalfNormal","Normal","Square-periodic",...
%     "tLocationScale","Uniform","Uniform-Mix","Weibull","Chisquare",...
%     "InverseGaussian","Trimodal-Normal","Stable",...
%     "Stable2","Stable3","Stable1","BirnbaumSaunders-Stable"];
% % % distributionVector = ["Weibull","Chisquare",...
% % %     "InverseGaussian","Trimodal-Normal","Stable",...
% % %     "Stable2","Stable3","Stable1","BirnbaumSaunders-Stable"];
% distributionVector = ["Gamma"];
% distributionVector = ["Generalized-Pareto","Square-periodic","Stable"];
% distributionVector = ["Stable"];
% distributionVector = ["Beta-a0p5-b0p5"];
% distributionVector = ["Stable","Stable1","Stable2","Stable3"];
distributionVector = ["Beta-a0p5-b0p5"];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Main Function Call Loop used to lable plot figures

% change code to multithreaded version
for j = 1:length(distributionVector)
    % Define plot vector for distributions from 0-1
    if distributionVector(j) == "Beta-a0p5-b1p5" ||...
            distributionVector(j) == "Beta-a2-b0p5" ||...
            distributionVector(j) == "Beta-a0p5-b0p5"
        lowLim = 0;
        upLim = 1;
        synx = linspace(lowLim,upLim,1000);
    else
        % Define plot vector for distribution from lowLim-upLim
        lowLim = 0;
        upLim = 10;
        synx = linspace(lowLim,upLim,1000);
    end
    % Current distribution name
    distributionName = distributionVector(j);
    % file name for actual distribution. "A_" puts at the top of the folder.
    fileNameAct = sprintf(['A_', char(distributionName),'_Act']);
    
    for i = 1:ntrials
        % Create vector of  samples
        sampleVec = samplesVector(minSamplesExp,maxSamplesExp,dataTypeflag,step);
        % initialize vector to track time of computation for sample size
        twallVec = [];
        tcpuVec = [];
        ticVec = [];
        sampleTrack = [];
        binData = [];
        T = [];
        
        for k = 1:length(sampleVec)
            interpEst = [];
            Ns = sampleVec(k);
            realx = linspace(lowLim,upLim,Ns);
            % p-vector definition for Rtree
            p = [1,0.5,1,0.33,2,ceil(0.0625*Ns^0.5),40];
            
            % Create fileName for each distribtuion
            fileName = sprintf(['D_', char(distributionName),'_T_','%d', '_S_','%d'],i, Ns);
            
            % Run DJ estimator and store data  ////////////////////////////
            if SEcallflag
                % file path name
                sendFileName = ['D_',char(distributionName),'\',char(fileName),'.txt'];
                
                %==========================================================
                % % % % % % % start of subsample estimate % % % % % % % % %
                %==========================================================
                tintial = cputime;
                
                
                % initial details for subsample
                sendFileName1 = ['D_',char(distributionName),'\',char(fileName),'.txt'];
                sample = importdata(sendFileName1);
                Ns = length(sample);
                % boot strap sample details
                bootNs = floor(percSample*Ns);
                numSample = numSubs;
                
                subSampPlot = [];
                Tsub = zeros(1,numSample);
                % loop to get estimate for bootstrapped samples
                for s = 1:numSample
                    disp('   ')
                    disp(['Sample type: ',char(fileName),'.txt'])
                    % Bootstrap
                    subSample = datasample(sample,bootNs,'Replace',false);
                    subSampPlot = [subSampPlot,subSample];
                    sendFileName = ['BS000',num2str(s),'.dat'];
                    dlmwrite(sendFileName,subSample,'Precision',12)
                    
                    [T0,DJ_x,DJ_pdf,DJ_cdf,DJ_u,DJ_SQR,nBlocks,Blacklist,Ns,binNs]...
                        = stitchPDF(fileName,sendFileName,savePNG,lowLim,upLim,p);
                    % store threshold value per subsample for given sample
                    % size
                    Tsub(1,s) = T0;
                    % store subsample estimate data
                    estimateData(s).x = DJ_x;
                    estimateData(s).pdf = DJ_pdf;
                    estimateData(s).u = DJ_u;
                    estimateData(s).sqr = DJ_SQR;
                end
                
                % AVERAGE SUBSAMPLES---------------------------------------
                % generate sample to interpolate with unique elements
                minX = zeros(size(estimateData,2),1);
                maxX = zeros(size(estimateData,2),1);
                for jj = 1:size(estimateData,2)
                    minX(jj) = min(estimateData(jj).x);
                    maxX(jj) = max(estimateData(jj).x);
                end
                % create interpolate x-range
                %                 xinterp = [max(minX):0.001:min(maxX)]';
                %^^^ not fine for heavey tailed distributions
                dxStep = 0.01;
                disp('---------------------------------------------------')
                disp(['num of data points from estimate: ',num2str(size(estimateData(1).x,2))])
                disp(['num of data points interpilated: ',num2str((min(maxX)-max(minX))/dxStep)])
                disp('---------------------------------------------------')
                % test that will interpilate ll estimates to the data
                % points of the first
                %                 xinterp = [max(minX):dxStep:min(maxX)]';
                %                 xinterp = linspace(max(minX),min(maxX),1000);
                %                 xinterp = linspace(max(minX),min(maxX),size(estimateData(1).x,2));
                xinterp = estimateData(1).x;
                % assign interpolated estimates into matrix columns
                for jj = 1:size(estimateData,2)
                    %
                    interpEst(:,jj) = interp1(estimateData(jj).x,estimateData(jj).pdf,xinterp);
                end
                % find average estimate
                AvgEst = mean(interpEst,2);
                % find average T per sample size and store subsamples
                T = [T,Tsub'];
                Tavg = mean(Tsub);
                % find z-score
                %                 if numSubs > 1
                %                     standDev = std(interpEst);
                %                     z = (interpEst - AvgEst)./standDev;
                %                     weight = 2*normcdf(z);
                %                     AvgPDF = mean(weight.*interpEst,2);
                %                     % compute estimate for estimated x-points
                %                     tempStruc.minVariance = 1;
                %                     for tt = 1:size(interpEst,1)
                %                         try
                %                             [failed, testX{tt}, testPDF{tt}, ~, ~,~] = EstimatePDF(z(tt,:),tempStruc);
                %
                %                             disp('NMEM failed')
                %                             disp(failed)
                %                         catch
                %                             warning(['Problem using function.  Assigning a value of -1.']);
                %                             testPDF{tt} = -1*ones(1,length(interpEst(tt,:)));
                %                             testX{tt} = -1*ones(1,length(interpEst(tt,:)));
                %                         end
                %                     end
                %                 end
                %----------------------------------------------------------
                tcpu = cputime-tintial;
                %==========================================================
                % % % % % % % % end of bootstrap estimate % % % % % % % % %
                %==========================================================
                
                % store time of computation
                tcpuVec = [tcpuVec,tcpu];
                sampleTrack = [sampleTrack,Ns];
                % read in actual distribution
                actSample = importdata(['D_',char(distributionName),'\',char(fileNameAct),'.txt']);
                pdfCurve = actSample(:,2);
                [pdfDistDJ] = distributionsChoices(distributionName,DJ_x,fileName,"off",precision,Ns);
                % Create final answer file
                fileName = sprintf(['stitchPDF_D_', char(distributionName),'_T_','%d', '_S_','%d'],i, Ns);
                Sanswer = [DJ_x',DJ_pdf'];
                dlmwrite([fileName,'.txt'],Sanswer, 'delimiter',' ','precision',12)
            end
            
            % SE PLOTS-----------------------------------------------------
            if SEplotflag
                %                 if numSubs > 1
                %                     figure('Name','Zscore Distribution')
                %                     hold on
                %                     for tt = 1:size(AvgPDF,1)
                %                         plot(testX{tt},testPDF{tt})
                %                     end
                %                     title(['Z-score_',char(distributionName)],'Interpreter','latex')
                %                     if savePNG
                %                         pngfile = strcat('Zscore_D_',char(distributionName),'T_',int2str(i),'S_',int2str(sampleVec(k)),'.png');
                %                         saveas(gcf,pngfile)
                %                         figfile = strcat('Zscore_D_',char(distributionName),'T_',int2str(i),'S_',int2str(sampleVec(k)),'.fig');
                %                         saveas(gcf,figfile)
                %                     end
                %                 end
                %----------------------------------------------------------
                figure('Name','Bootstrap Estimates')
                hold on
                l = zeros(1,2);
                plot(synx,pdfCurve,'--k')
                l(1) = plot(estimateData(1).x,estimateData(1).pdf,'Color',[.8,.8,.8],'DisplayName','Estimates for sNs');
                for s = 2:size(estimateData,2)
                    plot(estimateData(s).x,estimateData(s).pdf,'Color',[.8,.8,.8])
                end
                plot(xinterp,AvgEst,'Color',[0,0,0],'DisplayName','Average');
                %                 for next = 1:size(estimateData.x,1)
                %                 plot(estimateData(next).x,0.1*next*ones(1,size(estimateData(next).x,2)),'.m')
                %                 end
                
                if max(DJ_pdf) > 2
                    ylim([0,6])
                else
                    ylim([0,1])
                end
                xlim([lowLim,upLim])
                ylabel('$PDF$','Interpreter','latex')
                xlabel('x','Interpreter','latex')
                title(['bNs: ',num2str(bootNs),' Ns: ', num2str(Ns),' blocks: ',num2str(nBlocks),' numSubs: ',num2str(numSubs),' percSample: ',num2str(percSample)],'Interpreter','latex')
                legend(l(1))
                if savePNG
                    pngfile = strcat('boot_PDF_D_',char(distributionName),'T_',int2str(i),'S_',int2str(sampleVec(k)),'_',int2str(percSample),'_',int2str(numSubs),'.png');
                    saveas(gcf,pngfile)
                    figfile = strcat('boot_PDF_D_',char(distributionName),'T_',int2str(i),'S_',int2str(sampleVec(k)),'_',int2str(percSample),'_',int2str(numSubs),'.fig');
                    saveas(gcf,figfile)
                end
                
                % SQR plot function
                [boot_u,boot_sqr] = SQR(xinterp,AvgEst,sample);
                
                
                %----------------------------------------------------------
                figure('Name','Stitch Estimator: SQR')
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
                ldg2(1) = plot(boot_u,boot_sqr,'-k','DisplayName','BSSE');
                xlim([0,1])
                ylabel('$SQR$','Interpreter','latex')
                xlabel('u','Interpreter','latex')
                title(['SQR vs u','for Ns: ',int2str(Ns), ' and ',char(distributionName),' percSample: ',num2str(percSample),' numSubs: ',num2str(numSubs)],'Interpreter','latex');
                if savePNG
                    pngfile = strcat('SQR_D_',char(distributionName),'T_',int2str(ntrials),'S_',int2str(sampleVec(k)),'_',int2str(percSample),'_',int2str(numSubs),'.png');
                    saveas(gcf,pngfile)
                    figfile = strcat('SQR_D_',char(distributionName),'T_',int2str(ntrials),'S_',int2str(sampleVec(k)),'_',int2str(percSample),'_',int2str(numSubs),'.fig');
                    saveas(gcf,figfile)
                end
            end
        end
        
        size(T)
        size(sampleVec)
        
        figure('Name','Threshold per sample')
        hold on
        plot(sampleVec,Tavg,'-r')
        for i = 1:size(T,1)
            plot(sampleVec,T(i,:))
        end
        title('Track threshold per sample size')
        
        % sync computation times into matrix
        CtcpuT = horzcat(sampleTrack',tcpuVec')
        fileName = sprintf(['D_', char(distributionName),'_T_','%d'],i);
        dlmwrite(['CtcpuT',fileName,'_',int2str(percSample),'_',int2str(numSubs),'.txt'],CtcpuT,'delimiter',' ')
    end
end
toc
% end
% ^^ when function header is included