function [T0,sx,sPDF,sCDF,u,sqr,nBlocks,Blacklist,Ns,binNs] = stitchPDF(filename,sendFileName,savePNG,lowLim,upLim,p)

% stitchPDF switching board %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plotQQandSQR =              false; %<- true/false plot on/off QQ&SQR from stitchpdf
saveFIG =                   false; %<- true/false save plot figures on/off
% Script Detailed output %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
show_block_layout =         true; % true/false see details of block chain
visualize_interpolation =   false; % true/false  checking & debuging on/off
blockPdfs =                 true; % true/false plot on/off pdf for each block

%%                                                              user input      section 1
targetCoverage = 65;
targetCoverage = round(targetCoverage);
minNs = 10;                                    % minimum number of samples
maxNs = 100000000;                        % larger makes program very slow
%%                                                     some error checking      section 2
targetCoverage = sort(targetCoverage);
nTargets = length(targetCoverage);
if( nTargets < 1 )
    error('no targets specified');
end

% --------------------------------------------- universal scoring function
logLikelihood = LL();
xScore0 = logLikelihood(:,1);
yCoverage0 = logLikelihood(:,2);
xScore = xScore0;
yCoverage = 100*yCoverage0;

%%                                                          define weights      section 3
tL = targetCoverage - 2.5;
tH = targetCoverage + 2.5;
if( tL < 0.099999 )
    error('A target coverage value is too low: Lowest = 2.6');
end
if( tH > 95.000001 )
    error('A target coverage value is too high: Highest = 92.5');
end
% --------------------------------------------------------- assign weights
indx = 15:length(xScore)-15;
dx = xScore(indx+14) - xScore(indx-14);
dy = yCoverage(indx+14) - yCoverage(indx-14);
wScore = dy./dx;
xScore = xScore(15:end-15);
%plot(xScore,wScore);                                        % for checking
wt = size(targetCoverage);
for k=1:nTargets
    a = yCoverage - targetCoverage(k);
    [~,indx] = min( abs(a) );
    wt(k) = wScore(indx);
end
wtsum = sum(wt);
wt = wt/wtsum;

% figure('Name','Target Coverage')
% plot(targetCoverage,wt);                                   % for checking
% xlabel('Target Coverage')
% ylabel('Weight')
% %
% figure('Name','Likelihood')
% hold on
% % plot(xScore0,yCoverage0,'-r');                           % for checking
% plot(xScore,wScore,'-b');
% xlabel('xScore0')
% ylabel('yCoverage0')

%%                                          define qualitative descriptors      section 4
%---------------------------------------------------------- read in sample
FileName4Sample = sendFileName;
sample = sort(importdata(FileName4Sample));
x = sort(sample); 
Ns = length(sample);
binNs = 15;
if( Ns < minNs )
    error(['Sample size is too small. minNs = ',num2str(minNs)]);
elseif( Ns > maxNs )
    error(['Sample size is too large. maxNs = ',num2str(maxNs)]);
end

msgModelType = 'unknown distribution';
prefix = 'RD';

if( maxNs > 999999999999 )
    error('format error on prefix. do not let maxNs > 9999999');
else
    prefix = [prefix,num2str(Ns,'%012d')];
end

%--------------------------------------------------- display useful output
disp(['Sample file name: ',FileName4Sample]);
disp(['max sample point: ',num2str(max(sample))]);
disp(['min sample point: ',num2str(min(sample))]);
disp(['number of samples: Ns = ',num2str(Ns)]);
disp(['minimum pts per bin: binNs = ',num2str(binNs)]);
%%                          OBSERVATION WINDOW: CDF Boundary Probabilities      section 5
pL = 0.5/Ns;   % probability for data to be  left of window
pR = 0.5/Ns;   % probability for data to be right of window
pNorm = 1 - pL - pR;          % probability for data to fall within window
%%                                      partition data into primary blocks      section 6

[T0,j,nBlocks,kBlockLower,kBlockUpper,kList] = binWidthSize(Ns,binNs,sample,filename,savePNG,lowLim,upLim,p);

%%                                     create staggered (secondary) blocks      section 7
if( nBlocks > 2 )
    if( nBlocks == 3 )          % block size of three is a special case
        kBlockLower(1) = 1;
        kBlockUpper(1) = kList(1);
        kBlockLower(2) = round( ( 1  + kList(1) )/2 );
        kBlockUpper(2) = round( ( kList(1) + Ns )/2 );
        kBlockLower(3) = kList(1);
        kBlockUpper(3) = Ns;
    else
        % set index for staggered blocks (sb)
        nB = 1;
        % set edge positions for first sb
        kBlockLower(nB) = 1;
        kBlockUpper(nB) = kList(1);
        nB = nB + 1;
        % set edge positions for second sb as average first blocks (fb)
        kBlockLower(nB) = round( ( 1 + kList(1) )/2 );
        kBlockUpper(nB) = round( ( kList(1) + kList(2) )/2 );
        
        % Loop to decide other sb edge positions
        for i=1:j-2
            % update indexer for kBlock indices
            nB = nB + 1;
            % set edge positions as first block for i-th sb
            kBlockLower(nB) = kList(i);
            kBlockUpper(nB) = kList(i+1);
            % update indexer for kBlock indices
            nB = nB + 1;
            % set edge positions as average of fb positions for i+1-th sb
            kBlockLower(nB) = round( ( kList(i)  +  kList(i+1) )/2 );
            kBlockUpper(nB) = round( ( kList(i+1) + kList(i+2) )/2 );
        end
        nB = nB + 1;
        kBlockLower(nB) = kList(j-1);
        kBlockUpper(nB) = kList(j);
        nB = nB + 1;
        kBlockLower(nB) = round( ( kList(j-1)  +  kList(j) )/2 );
        kBlockUpper(nB) = round( ( kList(j)  +  Ns )/2 );
        
        kBlockLower(nBlocks) = kBlockUpper(end-2);
        kBlockUpper(nBlocks) = Ns;
    end
else                                               % => block size is 1
    kBlockLower = zeros(1,nBlocks);
    kBlockUpper = zeros(1,nBlocks);
    kBlockLower(1) = 1;
    kBlockUpper(1) = Ns;
end

%%                                                   Determine block sizes      section 8
blockSize = kBlockUpper - kBlockLower + 1;
flag = 1;
if( blockSize(1) < binNs )
    flag = -1;
    k = binNs - blockSize(1);
    %                 disp('k')
    %                 disp(k)
    tt = blockSize - k;
    %                 disp('tt')
    %                 disp(tt)
    for j=2:nBlocks
        if( tt(j) > binNs )
            jRef = j;
            if( jRef == nBlocks )
                error('All blocks are too small! Sample is too hard');
            end
            break;
        end
    end
    kBlockLower(1) = 1;
    kBlockUpper(1) = kBlockUpper(1) + k;
    for j=2:jRef-1
        kBlockUpper(j) = kBlockUpper(j) + k;
        kBlockLower(j) = kBlockLower(j) + k;
    end
    kBlockLower(jRef) = kBlockLower(jRef) + k;
    kBlockLower(jRef+1) = kBlockLower(jRef+1) + k;
end
if( blockSize(nBlocks) < binNs )
    flag = -1;
    k = binNs - blockSize(nBlocks);
    tt = blockSize - k;
    for j=nBlocks-1:-1:1
        if( tt(j) > binNs )
            jRef = j;
            if( jRef == 1 )
                error('All blocks are too small! Sample is too hard');
            end
            break;
        end
    end
    kBlockUpper(nBlocks) = Ns;
    kBlockLower(nBlocks) = kBlockLower(nBlocks) - k;
    for j=nBlocks-1:-1:jRef+1
        kBlockUpper(j) = kBlockUpper(j) - k;
        kBlockLower(j) = kBlockLower(j) - k;
    end
    kBlockUpper(jRef) = kBlockUpper(jRef) - k;
    kBlockUpper(jRef-1) = kBlockUpper(jRef-1) - k;
end
if( flag < 0 )
    blockSize = kBlockUpper - kBlockLower + 1;
end

%%                                get length scale & shifts for each block      section 9
blockScale = zeros(1,nBlocks);     % => length of x as the span of a block
blockShift = zeros(1,nBlocks);     % shift in x to center the span about 0
for b=1:nBlocks
    kL = kBlockLower(b);
    kU = kBlockUpper(b);
    blockScale(b) = 0.5*( x(kU) - x(kL) );          % REMEMBER x is sorted!
    %blockScale(b) = 0.1;                              % for debugging only
    blockShift(b) = 0.5*( x(kU) + x(kL) );          % REMEMBER x is sorted!
    %blockShift(b) = 0;                                % for debugging only
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if b == 9
        blockShift(b) = 0.5*( x(kU) + x(kL) );
    elseif b == 10
        blockShift(b) = 0.5*( x(kU) + x(kL) );
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
% % disp('BlockShift')
% % disp(blockShift)
% % pause
% ----------------------------------------------------------- for checking
if( show_block_layout )
    t = 1:nBlocks;
    figure('Name','BlockScale')
    plot(t,blockScale,'-k');
    title('Length of each block with respect to x variable');
    ylabel('Length of Block','Interpreter','latex')
    xlabel('Block Index','Interpreter','latex')
    title('Length Block with Respect to $x$','Interpreter','latex')
    if savePNG
        pngfile = strcat('BlockLength_',char(filename),'.png');
        saveas(gcf,pngfile)
        figfile = strcat('BlockLength_',char(filename),'.fig');
        saveas(gcf,figfile)
    end
    figure('Name','BlockSize')
    plot(t,blockSize,'-k');
    ylabel('Number of Data Points','Interpreter','latex')
    xlabel('Block Index','Interpreter','latex')
    title('Data Points per Block','Interpreter','latex')
    if savePNG
        pngfile = strcat('BlockSize_',char(filename),'.png');
        saveas(gcf,pngfile)
        figfile = strcat('BlockSize_',char(filename),'.fig');
        saveas(gcf,figfile)
    end
end
aveBlockSize = mean(blockSize);
maxBlockSize = max(blockSize);
minBlockSize = min(blockSize);
stdBlockSize = std(blockSize);
%%                                              report with error checking      section 10
if( nBlocks == 1 )
    disp('   number of blocks = 1');
    kL = 1;
    kR = Ns;
    disp(['   block 1: (kL,kR) = (',num2str(kL),',',num2str(kR),')']);
else
    disp(['  number of blocks = ',num2str(nBlocks)]);
    disp(['average block size = ',num2str(aveBlockSize)]);
    disp(['maximum block size = ',num2str(maxBlockSize)]);
    disp(['minimum block size = ',num2str(minBlockSize)]);
    disp(['std dev block size = ',num2str(stdBlockSize)]);
    %-------------------------------------------------------- for debugging
    if false
        % ----------------------------- write out details on 1st and 2nd blocks
        kL = kBlockLower(1);
        kR = kBlockUpper(1);
        disp(['   block 1: (kL,kR) = (',num2str(kL),',',num2str(kR),')']);
        % 2nd block
        kL = kBlockLower(2);
        kR = kBlockUpper(2);
        disp(['   block 2: (kL,kR) = (',num2str(kL),',',num2str(kR),')']);
        % ---------------------- write out details on 2nd to last & last blocks
        kL = kBlockLower(end-1);
        kR = kBlockUpper(end-1);
        disp(['block Ns-1: (kL,kR) = (',num2str(kL),',',num2str(kR),')']);
        kL = kBlockLower(end);
        kR = kBlockUpper(end);
        disp(['  block Ns: (kL,kR) = (',num2str(kL),',',num2str(kR),')']);
    else
        for index = 1:length(kBlockLower)
            kL = kBlockLower(index);
            kR = kBlockUpper(index);
            disp(['   block ',num2str(index),': (kL,kR) = (',num2str(kL),',',num2str(kR),')']);
            if index > 9
                space = '             ';
            elseif index > 99
                space = '              ';
            else
                space = '            ';
            end
            disp([space,'(xL,xR) = (',num2str(x(kL)),',',num2str(sample(kR)),')'])
            
        end
    end
end
disp(['  target coverage = ',num2str(targetCoverage)]);
%%                                              set up probability factors      section 11
nBs0 = blockSize - 1;
nBs0(1) = nBs0(1) + 0.5;
nBs0(nBlocks) = nBs0(nBlocks) + 0.5;
%%                                              divide problem into blocks      section 12
totalCounts = nBlocks*nTargets;
% --------------------------------------------------- commit to full Monte
% figure('Name','Test')
% hold on 
% for b=1:nBlocks
%     j1 = kBlockLower(b);
%     j2 = kBlockUpper(b);
%     disp(['x(',num2str(b),')'])
%     %     disp(size(x(j1+1:j2-1)))
%     
%     vizX = [x(j1),x(j2)];
%     vizY = b*[1,1];
%     plot(vizX,vizY,'-')
%     xlim([min(x),max(x)])
%     
%     blockCheckPrior{1,b}.X(1,:) = vizX;
%     blockCheckPrior{1,b}.Y(1,:) = vizY;
%     
% end
% 
% pause
for b=1:nBlocks
    j1 = kBlockLower(b);
    j2 = kBlockUpper(b);
    s = 1/blockScale(b);
    
    sample = s*( x(j1+1:j2-1) - blockShift(b) )';       % Note: x is sorted
    inFileName = [prefix,'b',num2str(b,'%03d'),'.dat'];
    dlmwrite(inFileName,sample,'delimiter','\n','precision',12);
end
%%                                                        run PDFestimator      section 13
tic
%         parameters.SURDtarget = targetCoverage(1);
%         parameters.SURDmin = 1;
%         parameters.SURDmax  = 100;
%         parameters.LagrangeMin = 1;
%         parameters.LagrangeMax = 200;
%         parameters.integrationPoints  = 300;
%         parameters.debug = true;
%         tempStruc.minVariance = 1;

parfor b=1:nBlocks
    baseName = [prefix,'b',num2str(b,'%03d')];
    inFileName = [baseName,'.dat'];
    MatPDFsample = importdata(inFileName);
    for t=1:nTargets
        tempStruc = struct('SURDtarget',targetCoverage(t));
        if b == 1
            tempStruc.highBound = 1;
            %             tempStruc.lowBound = [];
%             tempStruc.integrationPoints = 300;
        elseif b == nBlocks
            tempStruc.lowBound = -1;
            %             tempStruc.highBound = [];
%             tempStruc.integrationPoints = 2000;
        else
            tempStruc.lowBound = -1;
            tempStruc.highBound = 1;
%             tempStruc.integrationPoints = 300;
        end
        
%         tempStruc.debug = true;
%         tempStruc.SURDmin = 20;
%         tempStruc.SURDmax  = 80;
%         tempStruc.integrationPoints = 3000;
%         tempStruc.integrationPoints = 400;
        tempStruc.minVariance = 1;
        try
            [failed, targetBlock{t,b}.data(:,1), targetBlock{t,b}.data(:,2), targetBlock{t,b}.data(:,3), ~,~] = EstimatePDF(MatPDFsample,tempStruc);
        catch
            warning(['Problem using function.  Assigning a value of 0.',' t: ',num2str(t),' b: ',num2str(b)]);
            targetBlock{t,b}.data(:,1) = linspace(-1,1,100);
            targetBlock{t,b}.data(:,2) = 0*ones(100,1);
            targetBlock{t,b}.data(:,3) = 0*ones(100,1);
        end
    end
end

% % figure('Name','Test')
% % hold on 
% % for b=1:nBlocks
% %     j1 = kBlockLower(b);
% %     j2 = kBlockUpper(b);
% %     %     disp(['x(',num2str(b),')'])
% %     %     disp(size(x(j1+1:j2-1)))
% %     
% %     vizX1 = [x(j1),x(j2)];
% %     vizY1 = b*[1,1];
% %     
% %     sInv = blockScale(b);             % Note:  sInv = 1/s = 1/blockScale(b)
% %     s = 1/sInv;
% %     
% %     vizY2 = (b+0.2)*[1,1];
% %     vizX2 = [min(targetBlock{1,b}.data(:,1)),max(targetBlock{1,b}.data(:,1))];
% %     
% %     vizY2 = (nBlocks-b+0.2)*[1,1];
% %     vizX2 = sInv*vizX2 + blockScale(b);
% %     
% %     plot(vizX1,vizY1,'-r')
% %     plot(vizX2,vizY2,'-b')
% %     xlim([min(x),max(x)])
% % end
% pause

%%                  read files and combine different target PDFs per block      section 14
disp('   ');
disp('Combining PDFs for same block from different targets ... ');
disp('   ');
disp('Start----------------------------------------------Finish');
blockPDF = cell(1,nBlocks);
blockCDF = cell(1,nBlocks);
blockX = cell(1,nBlocks);
fileW = cell(1,totalCounts);
nW = 0;
b = 1;
updateIndex = 1;
indexList = [];
Blacklist = [];
while( b <= nBlocks )
    baseName = [prefix,'b',num2str(b,'%03d')];           % fixed per block
    
    % ----------------------------------------------------- get block data
    blockX{b} = targetBlock{1,b}.data(:,1);        % same x for all targets
    blockPDF{b} = wt(1)*targetBlock{1,b}.data(:,2);
    blockCDF{b} = wt(1)*targetBlock{1,b}.data(:,3);
    
    for t=2:nTargets
        blockPDF{b} = blockPDF{b} + wt(t)*targetBlock{t,b}.data(:,2);
        blockCDF{b} = blockCDF{b} + wt(t)*targetBlock{t,b}.data(:,3);
    end
    
    % --------------------------------------------------- scale block data
    sInv = blockScale(b);             % Note:  sInv = 1/s = 1/blockScale(b)
    s = 1/sInv;
    sWb = s*(nBs0(b)/Ns);  % Wb = nBs0/Ns = weight based on frequency count
    % ------------------------------------------------------ apply scaling
    blockPDF{b} = sWb*blockPDF{b}; % s on PDF is required to undo sInv on x
    blockX{b} = sInv*blockX{b} ...           % put x back to original units
        + blockShift(b); % translates x-original back to its location
    %indexList is used to evaluate goodblocks later in script
    indexList = [indexList,b];
    updateIndex = updateIndex + 1;
    
%     disp(['min(blockX{',num2str(b),'}): ', num2str(min(blockX{b}))])
%     disp(['max(blockX{',num2str(b),'}):', num2str(max(blockX{b}))])
    
    b = b + 1;
    
    
end
nBlocks = b-1;
fprintf('|<\n');
%%                                             report errors/discrepancies       section 15
if( nW > 0 )
    disp('  ');
    disp(['WARNING: ',num2str(nW),' files found with failed solutions']);
    disp('---------------------------------------- List of discrepancies');
    for d=1:nW
        disp( fileW{d} );
    end
    disp('--------------------------------------------------------------');
    disp('Note that due to these substitutions the weighting is incorrect');
    disp('  ');
    disp('  ');
end
%\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
% ========================================== Begin: stitch blocks together
%/////////////////////////////////////////////////////////////////////////
% ------------------------------ for debugging: quick look at the sections
if blockPdfs
    figure('Name','blockPdfs')
    hold on
    for b=1:length(indexList)
        plot( blockX{indexList(b)} , blockPDF{indexList(b)} )
    end
    plot(x(kBlockUpper),0.2*ones(length(kBlockUpper),1),'.b')
    plot(x(kBlockLower),0.1*ones(length(kBlockLower),1),'.r')
    ylabel('$\hat{f}(x)$','Interpreter','latex')
    xlabel('$x$','Interpreter','latex')
    title('$\hat{f}(x)$ per block','Interpreter','latex')
    if max( blockX{indexList(length(indexList))}) < 1.1
        xlim([0,1])
        ylim([0,6])
    else
        xlim([0,10])
        ylim([0,1])
    end
    if savePNG
        pngfile = strcat('pdfBlocks_',char(filename),'.png');
        saveas(gcf,pngfile)
        figfile = strcat('pdfBlocks_',char(filename),'.fig');
        saveas(gcf,figfile)
    end
end
%%                                     sync block data into a single range      section 16
sx = [];
sPDF = [];
for b=1:length(indexList)
    sx = [sx, blockX{indexList(b)}'];
    sPDF = [sPDF, blockPDF{indexList(b)}'];        % multiply blockPDF by s ?????
end
[sx,indx] = unique(sx);
sPDF = sPDF(indx);
%%                                              stitch the blocks together      section 17
%     [sPDF] = blockSticther(nBlocks,blockX,visualize_interpolation,sx,blockPDF,blockCDF,filename);
%     pause
% implement a lever-arm weigthing between two blocks across common overlap
% region. Use blockCDF as the weighting factors.
% CDFL1 = max(CDFL)   CDFL0 = min(CDFL)
% CDFR1 = max(CDFR)   CDFR0 = min(CDFR)
% u = (CDFL - CDFL0)/(CDFL1 - CDFL0)    v = (CDFR - CDFR0)/(CDFR1 - CDFR0)
% a = (1 - u)^2       b = v^2
% fL = a/(a + b)
% fR = b/(a + b)                    note that fL + fR = 1
%                                   note that fL --> 0 at most right point
%                                   note that fR --> 0 at most left  point
% for b=1:length(indexList)-1
for b=1:nBlocks-1
    
    xmin = min(blockX{indexList(b+1)});
    xmax = max(blockX{indexList(b)});
%     xmin = min(blockX{b+1});
%     xmax = max(blockX{b});
    
    Lnot = or( (sx < xmin) , (sx > xmax) );  % => outside of overlap region
    xStitch = sx(~Lnot);                          % within overlap region
    k0 = length(xStitch);
    
    if( k0 < 1 ) %51 <------------------------------ currently very small!
        disp(['    left block # = ',num2str(b)]);
        disp(['   right block # = ',num2str(b+1)]);
        disp(['            xmax = ',num2str(xmax)]);
        disp(['            xmin = ',num2str(xmin)]);
        disp(['         overlap = ',num2str(k0)]);
        error('overlap is too small!');
    end
    %----------------------------------------------------------------------
    [xb,indx] = unique(blockX{indexList(b)});
    yb = blockPDF{indexList(b)}(indx);
    PDFlower = interp1(xb,yb,xStitch);
    zb = blockCDF{indexList(b)}(indx);
    ub = interp1(xb,zb,xStitch);
    %----------------------------------------------------------------------
    [xb1,indx] = unique(blockX{indexList(b+1)});
    yb1 = blockPDF{indexList(b+1)}(indx);
    PDFupper = interp1(xb1,yb1,xStitch);
    zb1 = blockCDF{indexList(b+1)}(indx);
    vb = interp1(xb1,zb1,xStitch);
    %----------------------------------------------------------------------
    u0 = min(ub);
    u1 = max(ub);
    v0 = min(vb);
    v1 = max(vb);
    if v1 == v0
        v0 = v0 + 0.000001;
        % CHEAP FIX ^^^ <----------------------------------------------------------- (*)
    end
    u = (ub - u0)/(u1 - u0);
    % Some times v1=v0 which leads to v = nan
    v = (vb - v0)/(v1 - v0);
    power = 2;
    aLower = (1 - u).^(power);
    aUpper = v.^(power);
    temp = aLower + aUpper;
    f_lower = aLower./temp;
    f_upper = aUpper./temp;
    stitchPDF = PDFlower.*f_lower + PDFupper.*f_upper;
    % -------------------------------------------------------- error check
    if( visualize_interpolation )
        mb = length( blockX{indexList(b)} );
        mb1 = length( blockX{indexList(b+1)} );
        disp(['overlaping blocks ',num2str(b),' and ',num2str(b+1)]);
        disp([' left block length = ',num2str(mb)]);
        disp(['right block length = ',num2str(mb1)]);
        disp(['              xmin = ',num2str(xmin)]);
        disp(['              xmax = ',num2str(xmax)]);
        disp(['     overlap count = ',num2str(k0)]);
        disp('-------------------------------------------------------');
        disp('    ');
        %----------------------------------------------
        figure('Name','CDF & PDF Stitching')
        subplot(1,3,1)
        plot(xStitch,f_lower,'r');
        hold on;
        plot(xStitch,f_upper,'b');
        ylabel('$\hat{F}(x)$','Interpreter','latex')
        xlabel('x','Interpreter','latex')
        title('(a)','Interpreter','latex')
        %----------------------------------------------
        subplot(1,3,2)
        plot(xStitch,PDFlower,'r');
        hold on;
        plot(xStitch,PDFupper,'b');
        plot(xStitch,stitchPDF,'k');
        ylabel('$\hat{f}(x)$','Interpreter','latex')
        xlabel('x','Interpreter','latex')
        legend('$\hat{f}_k(x)$','$\hat{f}_{k+1}(x)$','CDF-Stitch $\hat{f}_s(x)$','Interpreter','latex')
        title('(b)','Interpreter','latex')
        %----------------------------------------------
        subplot(1,3,3)
        plot(xStitch,stitchPDF,'k');
        ylabel('$\hat{f}(x)$','Interpreter','latex')
        xlabel('x','Interpreter','latex')
        legend('$\hat{f}_k(x)$','$\hat{f}_{k+1}(x)$','CDF-Stitch $\hat{f}_s(x)$','Interpreter','latex')
        title('(c)','Interpreter','latex')
        if savePNG
            pngfile = strcat('CDF_Stitch_Index_',int2str(b),'_',char(filename),'.png');
            saveas(gcf,pngfile)
            figfile = strcat('CDF_Stitch_Index_',int2str(b),'_',char(filename),'.fig');
            saveas(gcf,figfile)
        end
    end
    %-------------------------------------------------- substitute into sx
    [~,indx] = min( abs(sx-xmin) );
    dk = k0 - 1;
    sPDF(indx:indx+dk) = stitchPDF;
    
end
%%                                                      normalize sPDF      section 18
% ^Technical note:  ^sPDF is already normalized but, we can force the
%                    boundaries as they should be, which is not automatic
% integrate assuming linear interpolation only
sCDF = zeros( size(sPDF) );
%----------------------------------------------------------------------
sCDF(1) = 0;
kmax = length(sCDF);
disp(['length(sCDF): ',num2str(kmax)])
for k=2:kmax
    fave = 0.5*( sPDF(k) + sPDF(k-1) );
    area = fave*( sx(k) - sx(k-1) );
    sCDF(k) = sCDF(k-1) + area;
end
temp = sCDF(kmax);
sCDF = pNorm*(sCDF/temp) + pL;  % recalling what pL, pR and pNorm are
%\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
% ============================================ END: stitch blocks together
%/////////////////////////////////////////////////////////////////////////

%%                                           calculate information to plot      section 19
sample = x(2:end-1);
% calculate & plot SQR  
%----------------------------------------------------------- adjust u range
sampleUpLim = max(sx);
sampleLoLim = min(sx);
[row, ~] = find(sample <= sampleUpLim & sample >= sampleLoLim);
u = interp1(sx,sCDF,sample(row));    % get corresponding u for each x in sample
% % % figure('Name','interp1 figure')
% % % hold on
% % % plot(sx,sCDF,'.r')
% % % plot(sample,1,'xb')
% % % plot(sample(row),u,'om')
% % % legend('sCDF','sample','u')
% % % pause
%--------------------------------------------------------------------------
%uref = (1:Ns-2)/(Ns - 1);              % both end points have been removed
uref = (1:size(sample(row),1))/(size(sample(row),1) - 1);
if( size(uref,1) ~= size(u,1) )
    u = u';
end
% --------------------------------------------------- get scaled residual
sqr = sqrt(Ns)*(u - uref); % normal formula has sqrt(Ns+2) but Ns -> Ns-2
% plot results
StichResultsPlot(plotQQandSQR,uref,u,msgModelType,Ns,prefix,saveFIG,sqr)
% --------------------------------------------------------------- finished
end