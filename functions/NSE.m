classdef NSE
    properties (Constant)
        max_bs = 100000;
        p1 = 1;
        p2 = 0.55;
        p3 = 1;
        p4 = 0.33;
        p5 = 2;
        p6 = 0.0625;
        p7 = 0.5;
        p8 = 40;
        p = [NSE.p1,NSE.p2,NSE.p3,NSE.p4,NSE.p5,NSE.p6,NSE.p7,NSE.p8];
        % boot strap parameters to be added????
    end
    properties
        % possible outputs
        failed;
        sx;
        sPDF;
        sCDF;
        u;
        sqr;
        nBlocks;
        N;
        binN;
        LG_max;
        LG_sum;
        T;
        BRlevel;
        BR0;
    end
    methods

%         function [obj.failed, obj.sx, obj.sPDF, obj.sCDF, u, sqr, obj.nBlocks, N, ...
%                 obj.binN, LG_max, LG_sum, T, obj.BRlevel, BR0]...
%                 = stitch(obj, inputSample) % for visualization

        function [obj] = stitch(obj, inputSample) % for visualization

            % initialze the obj.failed trip flag to false
            obj.failed = 0;

            % Script Detailed output %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if ~exist('plt_blocksize','var')
                plt_blocksize = false;
            end
            if ~exist('plt_blockpdf','var')
                plt_blockpdf = false;
            end
            if ~exist('plt_stitchpdf','var')
                plt_stitchpdf = false;
            end

            %%                                                              user input      section 1
            targetCoverage = 40;
            targetCoverage = round(targetCoverage);
            minNs = 10;                                    % minimum number of samples
            maxNs = 140000000;                        % larger makes program very slow
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

            % define weights      section 3
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
            wt = size(targetCoverage,1);
            for k=1:nTargets
                a = yCoverage - targetCoverage(k);
                [~,indx] = min( abs(a) );
                wt(k) = wScore(indx);
            end
            wtsum = sum(wt);
            wt = wt/wtsum;

            % define qualitative descriptors      section 4
            orginalSample = inputSample;

            sample = unique(orginalSample);
            numUnique = length(orginalSample) - length(sample);
            if numUnique > 0
                warning(['Number of dupliates: ', num2str(numUnique)])
            end

            x = sort(sample);
            obj.N = length(x);
            obj.binN = 15;
            if( obj.N < minNs )
                error(['Sample size is too small. minNs = ',num2str(minNs)]);
            elseif( obj.N > maxNs )
                error(['Sample size is too large. maxNs = ',num2str(maxNs)]);
            end

            if( maxNs > 999999999999 )
                error('Do not let maxNs > 9999999');
            end

            % OBSERVATION WINDOW: CDF Boundary Probabilities      section 5
            pL = 0.5/obj.N;   % probability for data to be  left of window
            pR = 0.5/obj.N;   % probability for data to be right of window
            pNorm = 1 - pL - pR;          % probability for data to fall within window

            % partition data into primary blocks      section 6
            % [j,obj.nBlocks,kBlockLower,kBlockUpper,kList,T,obj.BRlevel,BR0] = block_definition.bin_width_size(obj.N,obj.binN,x);

            obt_blocks = blocks;
            obt_blocks.sample = x;
            obt_blocks.binNs = obj.binN;
            [j,obj.nBlocks,kBlockLower,kBlockUpper,kList,obj.T,obj.BRlevel,obj.BR0] = obt_blocks.bin_width_size(obt_blocks);

            % create staggered (secondary) blocks      section 7
            if( obj.nBlocks > 2 )
                if( obj.nBlocks == 3 )          % block size of three is a special case
                    kBlockLower(1) = 1;
                    kBlockUpper(1) = kList(1);
                    kBlockLower(2) = round( ( 1  + kList(1) )/2 );
                    kBlockUpper(2) = round( ( kList(1) + obj.N )/2 );
                    kBlockLower(3) = kList(1);
                    kBlockUpper(3) = obj.N;
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
                    kBlockUpper(nB) = round( ( kList(j)  +  obj.N )/2 );

                    kBlockLower(obj.nBlocks) = kBlockUpper(end-2);
                    kBlockUpper(obj.nBlocks) = obj.N;
                end
            else                                               % => block size is 1
                kBlockLower = zeros(1,obj.nBlocks);
                kBlockUpper = zeros(1,obj.nBlocks);
                kBlockLower(1) = 1;
                kBlockUpper(1) = obj.N;
            end

            % Determine block sizes      section 8
            blockSize = kBlockUpper - kBlockLower + 1;
            flag = 1;
            if( blockSize(1) < obj.binN )
                flag = -1;
                k = obj.binN - blockSize(1);
                tt = blockSize - k;
                for j=2:obj.nBlocks
                    if( tt(j) > obj.binN )
                        jRef = j;
                        if( jRef == obj.nBlocks )
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
            if( blockSize(obj.nBlocks) < obj.binN )
                flag = -1;
                k = obj.binN - blockSize(obj.nBlocks);
                tt = blockSize - k;
                for j=obj.nBlocks-1:-1:1
                    if( tt(j) > obj.binN )
                        jRef = j;
                        if( jRef == 1 )
                            error('All blocks are too small! Sample is too hard');
                        end
                        break;
                    end
                end
                kBlockUpper(obj.nBlocks) = obj.N;
                kBlockLower(obj.nBlocks) = kBlockLower(obj.nBlocks) - k;
                for j=obj.nBlocks-1:-1:jRef+1
                    kBlockUpper(j) = kBlockUpper(j) - k;
                    kBlockLower(j) = kBlockLower(j) - k;
                end
                kBlockUpper(jRef) = kBlockUpper(jRef) - k;
                kBlockUpper(jRef-1) = kBlockUpper(jRef-1) - k;
            end
            if( flag < 0 )
                blockSize = kBlockUpper - kBlockLower + 1;
            end

            % get length scale & shifts for each block      section 9
            blockScale = zeros(1,obj.nBlocks);     % => length of x as the span of a block
            blockShift = zeros(1,obj.nBlocks);     % shift in x to center the span about 0
            for b=1:obj.nBlocks
                kL = kBlockLower(b);
                kU = kBlockUpper(b);
                blockScale(b) = 0.5*( x(kU) - x(kL) );          % REMEMBER x is sorted!
                blockShift(b) = 0.5*( x(kU) + x(kL) );          % REMEMBER x is sorted!
            end

            % ----------------------------------------------------------- for checking
            if( plt_blocksize )
                t = 1:obj.nBlocks;
                figure('Name','BlockScale')
                plot(t,blockScale,'-k');
                title('Length of each block with respect to x variable');
                ylabel('Length of Block','Interpreter','latex')
                xlabel('Block Index','Interpreter','latex')
                title('Length Block with Respect to $x$','Interpreter','latex')

                figure('Name','BlockSize')
                plot(t,blockSize,'-k');
                ylabel('Number of Data Points','Interpreter','latex')
                xlabel('Block Index','Interpreter','latex')
                title('Data Points per Block','Interpreter','latex')
            end

            % report with error checking      section 10
            if( obj.nBlocks == 1 )
                kL = 1;
                kR = obj.N;
            else
                for index = 1:length(kBlockLower)
                    kL = kBlockLower(index);
                    kR = kBlockUpper(index);
                end
            end

            % set up probability factors      section 11
            nBs0 = blockSize - 1;
            nBs0(1) = nBs0(1) + 0.5;
            nBs0(obj.nBlocks) = nBs0(obj.nBlocks) + 0.5;
            % divide problem into blocks      section 12
            totalCounts = obj.nBlocks*nTargets;
            % commit to full Monte
            for b=1:obj.nBlocks
                j1 = kBlockLower(b);
                j2 = kBlockUpper(b);
                s = 1/blockScale(b);

                sample = s*( x(j1+1:j2-1) - blockShift(b) )';       % Note: x is sorted
                orig_sample = x(j1+1:j2-1)';
                scaled_sample = x(j1+1:j2-1)';
                MatPDFsample{b} = sample;
            end
            % run PDFestimator      section 13
            tic

            % initialize vector to hold all lagrainge mutiplers per block
            LG = zeros(1,obj.nBlocks);
            parfor b=1:obj.nBlocks
                lagrange = [];
                for t=1:nTargets
                    tempStruc = struct('SURDtarget',targetCoverage(t));
                    try
                        [~, targetBlock{t,b}.data(:,1), targetBlock{t,b}.data(:,2), targetBlock{t,b}.data(:,3), ~,lagrange] = EstimatePDF(MatPDFsample{b});
                    catch
                        warning(['Problem using function.  Assigning a value of 0.',' t: ',num2str(t),' b: ',num2str(b)]);
                        lagrange = 0;
                        targetBlock{t,b}.data(:,1) = linspace(min(MatPDFsample{b}),max(MatPDFsample{b}),tempStruc.integrationPoints);
                        targetBlock{t,b}.data(:,2) = 0*ones(tempStruc.integrationPoints,1);
                        targetBlock{t,b}.data(:,3) = 0*ones(tempStruc.integrationPoints,1);
                    end
                    LG(1,b) = size(lagrange,1);
                end
            end
            obj.LG_max = max(LG);
            obj.LG_sum = sum(LG);

            % read files and combine different target PDFs per block      section 14
            blockPDF = cell(1,obj.nBlocks);
            blockCDF = cell(1,obj.nBlocks);
            blockX = cell(1,obj.nBlocks);

            b = 1;
            updateIndex = 1;
            indexList = [];
            while( b <= obj.nBlocks )
                % get block data
                blockX{b} = targetBlock{1,b}.data(:,1);        % same x for all targets
                blockPDF{b} = wt(1)*targetBlock{1,b}.data(:,2);
                blockCDF{b} = wt(1)*targetBlock{1,b}.data(:,3);

                for t=2:nTargets
                    blockPDF{b} = blockPDF{b} + wt(t)*targetBlock{t,b}.data(:,2);
                    blockCDF{b} = blockCDF{b} + wt(t)*targetBlock{t,b}.data(:,3);
                end

                % scale block data
                sInv = blockScale(b);             % Note:  sInv = 1/s = 1/blockScale(b)
                s = 1/sInv;
                sWb = s*(nBs0(b)/obj.N);  % Wb = nBs0/N = weight based on frequency count
                % apply scaling
                blockPDF{b} = sWb*blockPDF{b}; % s on PDF is required to undo sInv on x
                blockX{b} = sInv*blockX{b} ...           % put x back to original units
                    + blockShift(b); % translates x-original back to its location
                %indexList is used to evaluate goodblocks later in script
                indexList = [indexList,b];
                updateIndex = updateIndex + 1;

                b = b + 1;
            end

            obj.nBlocks = b-1;
            fprintf('|<\n');

            %\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
            % Begin: stitch blocks together
            %/////////////////////////////////////////////////////////////////////////
            % ------------------------------ for debugging: quick look at the sections
            if plt_blockpdf
                figure('Name','plt_blockpdf')
                hold on
                for b=1:length(indexList)
                    plot( blockX{indexList(b)} , blockPDF{indexList(b)} )
                end
                ylabel('$\hat{f}(x)$','Interpreter','latex')
                xlabel('$x$','Interpreter','latex')
                if max( blockX{indexList(length(indexList))}) < 1.1
                    ylim([0,6])
                else
                    ylim([0,1])
                end
            end

            % sync block data into a single range      section 16
            obj.sx = [];
            obj.sPDF = [];
            for b=1:length(indexList)
                obj.sx = [obj.sx, blockX{indexList(b)}'];
                obj.sPDF = [obj.sPDF, blockPDF{indexList(b)}'];        % multiply blockPDF by s ?????
            end
            [obj.sx,indx] = unique(obj.sx);
            obj.sPDF = obj.sPDF(indx);

            % stitch the blocks together      section 17

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
            for b=1:length(indexList)-1

                xmin = min(blockX{indexList(b+1)});
                xmax = max(blockX{indexList(b)});

                Lnot = or( (obj.sx < xmin) , (obj.sx > xmax) );  % => outside of overlap region
                xStitch = obj.sx(~Lnot);                          % within overlap region
                k0 = length(xStitch);
                if( k0 < 1 ) %51 <------------------------------ currently very small!
                    disp(['    left block # = ',num2str(b)]);
                    disp(['   right block # = ',num2str(b+1)]);
                    disp(['            xmax = ',num2str(xmax)]);
                    disp(['            xmin = ',num2str(xmin)]);
                    disp(['         overlap = ',num2str(k0)]);
                    warning('overlap is too small!');
                    obj.failed = 1;
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

                u = (ub - u0)/(u1 - u0);
                v = (vb - v0)/(v1 - v0);
                power = 2;
                aLower = (1 - u).^(power);
                aUpper = v.^(power);
                temp = aLower + aUpper;
                f_lower = aLower./temp;
                f_upper = aUpper./temp;
                stitchPDF = PDFlower.*f_lower + PDFupper.*f_upper;

                if plt_stitchpdf
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
                    %         legend('$\hat{f}_k(x)$','$\hat{f}_{k+1}(x)$','CDF-Stitch $\hat{f}_s(x)$',Interpreter='latex')
                    title('(b)','Interpreter','latex')
                    %----------------------------------------------
                    subplot(1,3,3)
                    plot(xStitch,stitchPDF,'k');
                    ylabel('$\hat{f}(x)$','Interpreter','latex')
                    xlabel('x','Interpreter','latex')
                    %         legend('$\hat{f}_k(x)$','$\hat{f}_{k+1}(x)$','CDF-Stitch $\hat{f}_s(x)$',Interpreter='latex')
                    title('(c)','Interpreter','latex')
                end


                [~,indx] = min( abs(obj.sx-xmin) );
                dk = k0 - 1;
                obj.sPDF(indx:indx+dk) = stitchPDF;

            end

            % normalize obj.sPDF      section 18
            % ^Technical note:  ^obj.sPDF is already normalized but, we can force the
            %                    boundaries as they should be, which is not automatic
            % integrate assuming linear interpolation only
            obj.sCDF = zeros( size(obj.sPDF) );
            obj.sCDF(1) = 0;
            kmax = length(obj.sCDF);
            for k=2:kmax
                fave = 0.5*( obj.sPDF(k) + obj.sPDF(k-1) );
                area = fave*( obj.sx(k) - obj.sx(k-1) );
                obj.sCDF(k) = obj.sCDF(k-1) + area;
            end
            temp = obj.sCDF(kmax);
            obj.sCDF = pNorm*(obj.sCDF/temp) + pL;  % recalling what pL, pR and pNorm are
            %\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
            % END: stitch blocks together
            %/////////////////////////////////////////////////////////////////////////

            % calculate information to plot      section 19
            sample = x(2:end-1);
            % calculate & plot SQR

            % adjust u range
            sampleUpLim = max(obj.sx);
            sampleLoLim = min(obj.sx);
            [row, ~] = find(sample <= sampleUpLim & sample >= sampleLoLim);
            obj.u = interp1(obj.sx,obj.sCDF,sample(row));    % get corresponding u for each x in sample
            uref = (1:size(sample(row),1))/(size(sample(row),1) - 1);
            if( size(uref,1) ~= size(obj.u,1) )
                obj.u = obj.u';
            end


            % get scaled residual
            obj.sqr = sqrt(obj.N)*(obj.u - uref); % normal formula has sqrt(N+2) but N -> N-2
            [obj.u,obj.sqr] = misc_functions.sqr(obj.sx,obj.sPDF,inputSample);
        end
    end
end