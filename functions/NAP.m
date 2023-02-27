classdef NAP
    properties (Constant)
        p1 = 1;
        p2 = 1/5;
        p3 = 1;
        p4 = 1/3;
        p5 = 2;
        p6 = 0.0625;
        p7 = 0.5;
        p8 = 40;
    end
    properties (SetAccess = public)
        max_bs=1e5;
    end
    properties
        % possible outputs
        p = [NAP.p1,NAP.p2,NAP.p3,NAP.p4,NAP.p5,NAP.p6,NAP.p7,NAP.p8];
        failed;
        sx;
        sPDF;
        sCDF;
        u;
        sqr;
        nBlocks;
        N;
        binN;
        T;
        BRlevel;
        BR0;
        % estimate data
        blocks_x;
        blocks_pdf;
        blocks_cdf;
        block_indx;
        % meta data
        LG_vals;
        LG;
        LG_max;
        LG_sum;
        block_size;
        block_scale;
    end
    methods
        function [obj] = stitch(obj, inputSample, serial)

            % initialze the obj.failed trip flag to false
            obj.failed = 0;

            % user input
            targetCoverage = round(70);

            % minimum number of allowed sample size
            minNs = 50; 

            % some error checking
            targetCoverage = sort(targetCoverage);
            nTargets = length(targetCoverage);

            % universal scoring function
            logLikelihood = utils.likelihood();
            xScore0 = logLikelihood(:,1);
            yCoverage0 = logLikelihood(:,2);
            xScore = xScore0;
            yCoverage = 100*yCoverage0;

            % define weights
            tL = targetCoverage - 2.5;
            tH = targetCoverage + 2.5;
            if( tL < 0.099999 )
                error('A target coverage value is too low: Lowest = 2.6');
            end
            if( tH > 95.000001 )
                error('A target coverage value is too high: Highest = 92.5');
            end

            % assign weights
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

            % define qualitative descriptors
            orginalSample = inputSample;

            sample = unique(orginalSample);
            numUnique = length(orginalSample) - length(sample);
            if numUnique > 0
                warning(['Number of dupliates: ', num2str(numUnique)])
            end

            % sort sample and delta-x array
            x = sort(sample);
            dxn = sample(2:end) - sample(1:end-1);
            dxns = sort(dxn);

            obj.N = length(x);
            if( obj.N < minNs )
                error(['Sample size is too small. minNs = ',num2str(minNs)]);
            end
            obj.binN = 15;

            % OBSERVATION WINDOW: CDF Boundary Probabilities
            % probability for data to be  left of window
            pL = 0.5/obj.N;
            % probability for data to be right of window
            pR = 0.5/obj.N;
            % probability for data to fall within window
            pNorm = 1 - pL - pR;

            % partition data into primary blocks
            obt_blocks = blocks;
            obt_blocks.sample = x;
            obt_blocks.dx = dxn;
            obt_blocks.dxs = dxns;
            obt_blocks.binNs = obj.binN;
            [j,obj.nBlocks,kBlockLower,kBlockUpper,kList,obj.T,...
                obj.BRlevel,obj.BR0] ...
                = obt_blocks.bin_width_size(obt_blocks);

            % create staggered (secondary) blocks
            if( obj.nBlocks > 2 )
                % block size of three is a special case
                if( obj.nBlocks == 3 )
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
                    % set edge positions for second sb as average first 
                    % blocks (fb)
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
            else % => block size is 1
                kBlockLower = zeros(1,obj.nBlocks);
                kBlockUpper = zeros(1,obj.nBlocks);
                kBlockLower(1) = 1;
                kBlockUpper(1) = obj.N;
            end

            % Determine block sizes
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

            % get length scale & shifts for each block
            blockScale = zeros(1,obj.nBlocks);     % => length of x as the span of a block
            blockShift = zeros(1,obj.nBlocks);     % shift in x to center the span about 0
            for b=1:obj.nBlocks
                kL = kBlockLower(b);
                kU = kBlockUpper(b);
                blockScale(b) = 0.5*( x(kU) - x(kL) );          % REMEMBER x is sorted!
                blockShift(b) = 0.5*( x(kU) + x(kL) );          % REMEMBER x is sorted!
            end

            % track meta data
            obj.block_size = blockSize;
            obj.block_scale = blockScale;

            % report with error checking
            if( obj.nBlocks == 1 )
                kL = 1;
                kR = obj.N;
            else
                for index = 1:length(kBlockLower)
                    kL = kBlockLower(index);
                    kR = kBlockUpper(index);
                end
            end

            % set up probability factors
            nBs0 = blockSize - 1;
            nBs0(1) = nBs0(1) + 0.5;
            nBs0(obj.nBlocks) = nBs0(obj.nBlocks) + 0.5;
            % divide problem into blocks
            totalCounts = obj.nBlocks*nTargets;


            % initialize tempStruc for parfor
            bounds = cell(1,obj.nBlocks);
            % commit to full Monte
            for b=1:obj.nBlocks
                j1 = kBlockLower(b);
                j2 = kBlockUpper(b);
                s = 1/blockScale(b);

                sample = s*( x(j1+1:j2-1) - blockShift(b) )';       % Note: x is sorted
                MatPDFsample{b} = sample;

                if b == 1
                    tempStruc.highBound = 1;
                elseif b == obj.nBlocks
                    tempStruc.lowBound = -1;
                else
                    tempStruc.lowBound = -1;
                    tempStruc.highBound = 1;
                end
                bounds{b} = tempStruc;
            end
            % run PDFestimator
            
            % initialize vector to hold all lagrainge mutiplers per block
            LG = zeros(1,obj.nBlocks);
            LG_vals = cell(1,obj.nBlocks);

            % conditionally execture parfor or for loop
            if serial
                parforArg = 0;
            else
                parforArg = Inf;
            end
            parfor (b=1:obj.nBlocks, parforArg)
                lagrange = [];
                for t=1:nTargets
                    try
                        [~, targetBlock{t,b}.data(:,1), targetBlock{t,b}.data(:,2), targetBlock{t,b}.data(:,3), ~,lagrange] = EstimatePDF(MatPDFsample{b}, bounds{b});    
                    catch
                        warning(['Problem using function.  Assigning a value of 0.',' t: ',num2str(t),' b: ',num2str(b)]);
                        lagrange = 0;
                        targetBlock{t,b}.data(:,1) = linspace(min(MatPDFsample{b}),max(MatPDFsample{b}),100);
                        targetBlock{t,b}.data(:,2) = 0*ones(100,1);
                        targetBlock{t,b}.data(:,3) = 0*ones(100,1);
                    end
                    LG(1,b) = size(lagrange,1);
                    LG_vals{b} = lagrange;
                end
            end

            % track meta data
            obj.LG_vals = LG_vals;
            obj.LG = LG;
            obj.LG_max = max(LG);
            obj.LG_sum = sum(LG);

            % read files and combine different target PDFs per block
            blockPDF = cell(1,obj.nBlocks);
            blockCDF = cell(1,obj.nBlocks);
            blockX = cell(1,obj.nBlocks);
            blockMat = cell(1,obj.nBlocks);

            b = 1;
            updateIndex = 1;
            indexList = [];
            while( b <= obj.nBlocks )
                % get block data
                blockX{b} = targetBlock{1,b}.data(:,1); % same x for all targets
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

            %\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
            % Begin: stitch blocks together
            %//////////////////////////////////////////////////////////////
            obj.blocks_x = blockX;
            obj.blocks_pdf = blockPDF;
            obj.blocks_cdf = blockCDF;
            obj.block_indx = indexList;

            % sync block data into a single range
            obj.sx = [];
            obj.sPDF = [];
            for b=1:length(indexList)
                obj.sx = [obj.sx, blockX{indexList(b)}'];
                obj.sPDF = [obj.sPDF, blockPDF{indexList(b)}'];
            end
            [obj.sx,indx] = unique(obj.sx);
            obj.sPDF = obj.sPDF(indx);

            % stitch the blocks together
            for b=1:length(indexList)-1

                xmin = min(blockX{indexList(b+1)});
                xmax = max(blockX{indexList(b)});
                xmin_list = mink(blockX{indexList(b+1)},5);
                xmax_list = maxk(blockX{indexList(b)},5);

                Lnot = or( (obj.sx < xmin) , (obj.sx > xmax) );  
                xStitch = obj.sx(~Lnot); % within overlap region
                k0 = length(xStitch);

                % overlap is too small
                if( k0 < 1 )
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

                 %----------------------------------------------------------------------
                [xb1,indx] = unique(blockX{indexList(b+1)});
                yb1 = blockPDF{indexList(b+1)}(indx);
                PDFupper = interp1(xb1,yb1,xStitch);

                power=2;
                u = 0: (1/(length(xStitch)-1)) : 1;
                aLower = (1 - u).^(power);
                aUpper = u.^(power);
                temp = aLower + aUpper;
                f_lower = aLower./temp;
                f_upper = aUpper./temp;
                stitchPDF = PDFlower.*f_lower + PDFupper.*f_upper;

                [~,indx] = min( abs(obj.sx-xmin) );
                dk = k0 - 1;
                obj.sPDF(indx:indx+dk) = stitchPDF;

            end

            % normalize obj.sPDF
            % Technical note:  
            % obj.sPDF is already normalized but, we can force the
            % boundaries as they should be, which is not automatic
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
            % recalling what pL, pR and pNorm are
            obj.sCDF = pNorm*(obj.sCDF/temp) + pL;
            %\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
            % END: stitch blocks together
            %//////////////////////////////////////////////////////////////

            % calculate SQR
            sample = x(2:end-1);

            % adjust u range
            sampleUpLim = max(obj.sx);
            sampleLoLim = min(obj.sx);
            [row, ~] = find(sample <= sampleUpLim & sample >= sampleLoLim);
            % get corresponding u for each x in sample
            obj.u = interp1(obj.sx,obj.sCDF,sample(row));
            uref = (1:size(sample(row),1))/(size(sample(row),1) - 1);
            if( size(uref,1) ~= size(obj.u,1) )
                obj.u = obj.u';
            end


            % get scaled residual
            % normal formula has sqrt(N+2) but N -> N-2
            obj.sqr = sqrt(obj.N)*(obj.u - uref);
            [obj.u,obj.sqr] = utils.sqr(obj.sx,obj.sPDF,inputSample);
        end
    end
end