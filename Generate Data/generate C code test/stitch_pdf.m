function [failed,sx,sPDF,sCDF,u,sqr,nBlocks,Blacklist,Ns,binNs, LG_max, LG_sum,T,BRlevel,BR0] = stitch_pdf(inputSample,filename,sendFileName,lowLim,upLim,p)

failed = 0;

% stitchPDF switching board %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plotQQandSQR =              false; %<- true/false plot on/off QQ&SQR from stitchpdf
saveFIG =                   true; %<- true/false save plot figures on/off
% Script Detailed output %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
show_block_layout =         false; % true/false see details of block chain
visualize_interpolation =   false; % true/false  checking & debuging on/off
blockPdfs =                 false; % true/false plot on/off pdf for each block

%%                                                              user input      section 1
targetCoverage = 40;%65;
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
logLikelihood = misc_functions.likelihood();
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
wt = size(targetCoverage,1);
for k=1:nTargets
    a = yCoverage - targetCoverage(k);
    [~,indx] = min( abs(a) );
    wt(k) = wScore(indx);
end
wtsum = sum(wt);
wt = wt/wtsum;

%%                                          define qualitative descriptors      section 4
%---------------------------------------------------------- read in sample
FileName4Sample = sendFileName;
orginalSample = inputSample;

% inFileName = ['fullSample_',filename,'.dat'];
% dlmwrite(inFileName,orginalSample,'delimiter','\n','precision',12);

sample = unique(orginalSample);
numUnique = length(orginalSample) - length(sample);

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
% disp(['Sample file name: ',FileName4Sample]);
% disp(['max sample point: ',num2str(max(sample))]);
% disp(['min sample point: ',num2str(min(sample))]);
% disp(['number of samples: Ns = ',num2str(Ns)]);
% disp(['minimum pts per bin: binNs = ',num2str(binNs)]);
%%                          OBSERVATION WINDOW: CDF Boundary Probabilities      section 5
pL = 0.5/Ns;   % probability for data to be  left of window
pR = 0.5/Ns;   % probability for data to be right of window
pNorm = 1 - pL - pR;          % probability for data to fall within window
%%                                      partition data into primary blocks      section 6

[j,nBlocks,kBlockLower,kBlockUpper,kList,T,BRlevel,BR0] = block_definition.bin_width_size(Ns,binNs,x,filename,lowLim,upLim,p);

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
    tt = blockSize - k;
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
    blockShift(b) = 0.5*( x(kU) + x(kL) );          % REMEMBER x is sorted!
end

%     inFileName = ['Kblocks_',filename,'.dat'];
%     dlmwrite(inFileName,[kBlockLower,kBlockUpper],'delimiter','\n','precision',12);


aveBlockSize = mean(blockSize);
maxBlockSize = max(blockSize);
minBlockSize = min(blockSize);
stdBlockSize = std(blockSize);
%%                                              report with error checking      section 10
if( nBlocks == 1 )
    kL = 1;
    kR = Ns;
else
    if false
        % ----------------------------- write out details on 1st and 2nd blocks
        kL = kBlockLower(1);
        kR = kBlockUpper(1);
        % 2nd block
        kL = kBlockLower(2);
        kR = kBlockUpper(2);
        kL = kBlockLower(end-1);
        kR = kBlockUpper(end-1);
        kL = kBlockLower(end);
        kR = kBlockUpper(end);
    else
        for index = 1:length(kBlockLower)
            kL = kBlockLower(index);
            kR = kBlockUpper(index);
            if index > 9
                space = '             ';
            elseif index > 99
                space = '              ';
            else
                space = '            ';
            end
        end
    end
end
% disp(['  target coverage = ',num2str(targetCoverage)]);
%%                                              set up probability factors      section 11
nBs0 = blockSize - 1;
nBs0(1) = nBs0(1) + 0.5;
nBs0(nBlocks) = nBs0(nBlocks) + 0.5;
%%                                              divide problem into blocks      section 12
totalCounts = nBlocks*nTargets;
% --------------------------------------------------- commit to full Monte

for b=1:nBlocks
    j1 = kBlockLower(b);
    j2 = kBlockUpper(b);
    s = 1/blockScale(b);
    
    sample = s*( x(j1+1:j2-1) - blockShift(b) )';       % Note: x is sorted
    orig_sample = x(j1+1:j2-1)';
    scaled_sample = x(j1+1:j2-1)';
    MatPDFsample{b} = sample;
    
end
%%                                                        run PDFestimator      section 13
tic


% initialize vector to hold all lagrainge mutiplers per block
LG = zeros(1,nBlocks);
parfor b=1:nBlocks
    lagrange = [];
    baseName = [prefix,'b',num2str(b,'%03d')];
    inFileName = [baseName,'.dat'];
    
    for t=1:nTargets
        tempStruc = struct('SURDtarget',targetCoverage(t));
        
        tempStruc.SURDtarget = targetCoverage(1);
        tempStruc.SURDmin = 5;
        tempStruc.SURDmax  = 100;
        tempStruc.LagrangeMin = 1;
        tempStruc.LagrangeMax = 200;
        tempStruc.integrationPoints  = 5000;
        tempStruc.integrationPoints  = 1000;
        tempStruc.scoreType = 'QZ';
        tempStruc.minVariance = false;
        tempStruc.adaptiveDx = false;
        
        [failed, targetBlock{t,b}.data(:,1), targetBlock{t,b}.data(:,2), targetBlock{t,b}.data(:,3), ~,lagrange] = EstimatePDF(MatPDFsample{b},tempStruc);
        
        LG(1,b) = size(lagrange,1);
        
    end
end

LG_max = max(LG);
LG_sum = sum(LG);

%%                  read files and combine different target PDFs per block      section 14

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
    
    % save block data for debugging
    %     dlmwrite([filename,'_blockX_',num2str(b),'.dat'],blockX{b},'Precision',12)
    %     dlmwrite([filename,'_blockPDF_',num2str(b),'.dat'],blockPDF{b},'Precision',12)
    %     dlmwrite([filename,'_blockCDF_',num2str(b),'.dat'],blockCDF{b},'Precision',12)
    
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

%%                                     sync block data into a single range      section 16
sx = [];
sPDF = [];
for b=1:length(indexList)
    sx = [sx, blockX{indexList(b)}'];
    sPDF = [sPDF, blockPDF{indexList(b)}'];        % multiply blockPDF by s ?????
end
[sx,indx] = unique(sx);
sPDF = sPDF(indx);

% disp('sPDF')
%  sum(~isfinite(sPDF))

%%                                              stitch the blocks together      section 17

for b=1:length(indexList)-1
    
    xmin = min(blockX{indexList(b+1)});
    xmax = max(blockX{indexList(b)});
    
    Lnot = or( (sx < xmin) , (sx > xmax) );  % => outside of overlap region
    xStitch = sx(~Lnot);                          % within overlap region
    k0 = length(xStitch);
    
    if( k0 < 1 ) %51 <------------------------------ currently very small!
        disp(['    left block # = ',num2str(b)]);
        disp(['   right block # = ',num2str(b+1)]);
        disp(['            xmax = ',num2str(xmax)]);
        disp(['            xmin = ',num2str(xmin)]);
        disp(['         overlap = ',num2str(k0)]);
        failed = 1;
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
        % CHEAP FIX for failed blocks^^^ <----------------------------------------------------------- (*)
    end
    if u1 == u0
        u0 = u0 + 0.000001;
        % CHEAP FIX for failed blocks ^^^ <----------------------------------------------------------- (*)
    end
    u = (ub - u0)/(u1 - u0);
    % Some times v1=v0 or u1=u0 which leads to division by zero
    v = (vb - v0)/(v1 - v0);
    power = 2;
    aLower = (1 - u).^(power);
    aUpper = v.^(power);
    temp = aLower + aUpper;
    f_lower = aLower./temp;
    f_upper = aUpper./temp;
    stitchPDF = PDFlower.*f_lower + PDFupper.*f_upper;
    % -------------------------------------------------------- error check
    
    %-------------------------------------------------- substitute into sx
    [~,indx] = min( abs(sx-xmin) );
    dk = k0 - 1;
    sPDF(indx:indx+dk) = stitchPDF;
    
end

if sum(~isfinite(sPDF)) || sum(~isfinite(sx))
    sum(~isfinite(sPDF))
    % pdf trigger non-finite values when estimate fails
    sum(~isfinite(sx))
    
end



%%                                                      normalize sPDF      section 18
% ^Technical note:  ^sPDF is already normalized but, we can force the
%                    boundaries as they should be, which is not automatic
% integrate assuming linear interpolation only
sCDF = zeros( size(sPDF) );

%----------------------------------------------------------------------
sCDF(1) = 0;
kmax = length(sCDF);
% disp(['length(sCDF): ',num2str(kmax)])
for k=2:kmax
    fave = 0.5*( sPDF(k) + sPDF(k-1) );
    area = fave*( sx(k) - sx(k-1) );
    sCDF(k) = sCDF(k-1) + area;
end
temp = sCDF(kmax);
sCDF = pNorm*(sCDF/temp) + pL;  % recalling what pL, pR and pNorm are

if sum(~isfinite(sCDF)) || sum(~isfinite(sx))
    sum(~isfinite(sCDF))
    % pdf trigger non-finite values when estimate fails
    sum(~isfinite(sx))
    %pause
end

sCDF(1) = 0;
kmax = length(sCDF);
% disp(['length(sCDF): ',num2str(kmax)])
for k=2:kmax
    fave = 0.5*( sPDF(k) + sPDF(k-1) );
    area = fave*( sx(k) - sx(k-1) );
    sCDF(k) = sCDF(k-1) + area;
end
temp = sCDF(kmax);
sCDF = pNorm*(sCDF/temp) + pL;

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
%--------------------------------------------------------------------------
%uref = (1:Ns-2)/(Ns - 1);              % both end points have been removed
uref = (1:size(sample(row),1))/(size(sample(row),1) - 1);
if( size(uref,1) ~= size(u,1) )
    u = u';
end


% --------------------------------------------------- get scaled residual
sqr = sqrt(Ns)*(u - uref); % normal formula has sqrt(Ns+2) but Ns -> Ns-2
% plot results

[u,sqr] = misc_functions.sqr(sx,sPDF,inputSample);
% --------------------------------------------------------------- finished
end