classdef NAP
    properties (Constant)
        p1 = 1;
        p2 = 0.55;
        p3 = 1;
        p4 = 0.33;
        p5 = 2;
        p6 = 0.0625;
        p7 = 0.5;
        p8 = 40;

        mod_estimate = true;
        trgt_SURD = 20;
        smooth = 100;
        max_lrg_mltply = 100;

    end
    properties (SetAccess = public)
        max_bs=1e5;
        minN = 500;
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

            orginalSample = inputSample;

            sample = unique(orginalSample);
            numUnique = length(orginalSample) - length(sample);
            if numUnique > 0
                warning(['Number of dupliates: ', num2str(numUnique)])
            end

            x = sort(sample);
            obj.N = length(x);

            if length(inputSample) >= obj.minN
                % initialze the obj.failed trip flag to false
                obj.failed = 0;

                % Script Detailed output %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                save_figs = false;
                if ~exist('plt_stitchpdf','var')
                    plt_stitchpdf = false;
                end
                if ~exist('plt_blocksize','var')
                    plt_blocksize = false;
                end
                if ~exist('plt_blockpdf','var')
                    plt_blockpdf = false;
                end
                if ~exist('plt_stitchpdf','var')
                    plt_stitchpdf = false;
                end

%                 plt_stitchpdf = true;

                %%                                                              user input      section 1
                targetCoverage = obj.trgt_SURD;
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
                logLikelihood = utils.likelihood();
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

                n_s = length(sample);
                dxn = sample(2:end) - sample(1:end-1);
                dxns = sort(dxn);
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

                obt_blocks = blocksJ;
                obt_blocks.sample = x;
                obt_blocks.dx = dxn;
                obt_blocks.dxs = dxns;
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

                % track meta data
                obj.block_size = blockSize;
                obj.block_scale = blockScale;

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


                % initialize tempStruc for parfor
                bounds = cell(1,obj.nBlocks);
                % commit to full Monte
                for b=1:obj.nBlocks
                    j1 = kBlockLower(b);
                    j2 = kBlockUpper(b);
                    s = 1/blockScale(b);
                    %
                    sample = s*( x(j1+1:j2-1) - blockShift(b) )';       % Note: x is sorted
                    MatPDFsample{b} = sample;
                    MatPDFsample_notscaled{b} = x(j1+1:j2-1);

                    if b == 1
                        tempStruc.highBound = 1;
                    elseif b == obj.nBlocks
                        tempStruc.lowBound = -1;
                    else
                        tempStruc.lowBound = -1;
                        tempStruc.highBound = 1;
                    end

                    tempStruc.SURDtarget = obj.trgt_SURD;
                    if obj.mod_estimate
                        tempStruc.smooth = obj.smooth;
                        tempStruc.LagrangeMax = obj.max_lrg_mltply;
                    end


%                     tempStruc.partition = 0;                   % {any number, zero for no partitioning}
                    tempStruc.outlierCutoff = 0;               % {number >= 0, zero to keep all outliers}

                    bounds{b} = tempStruc;
                end


                % run PDFestimator      section 13

                % initialize vector to hold all lagrainge mutiplers per block
                LG = zeros(1,obj.nBlocks);
                LG_vals = cell(1,obj.nBlocks);

                % conditionally execture parfor or for loop
                if serial
                    parforArg = 0;
                else
                    parforArg = Inf;
                end
%                 parfor (b=1:obj.nBlocks, parforArg)

                if serial
                    for b=1:obj.nBlocks
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
                else
                    ticBytes(gcp);
                    parfor b=1:obj.nBlocks
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
                    mem_bytes = tocBytes(gcp);

                    fprintf('Mega Bytes per parfor: %.2i MB\n', mem_bytes*1e-6)
                end

                % track meta data
                obj.LG_vals = LG_vals;
                obj.LG = LG;
                obj.LG_max = max(LG);
                obj.LG_sum = sum(LG);

                % read files and combine different target PDFs per block      section 14
                blockPDF = cell(1,obj.nBlocks);
                blockCDF = cell(1,obj.nBlocks);
                blockX = cell(1,obj.nBlocks);
                blockMat = cell(1,obj.nBlocks);

                b = 1;
                updateIndex = 1;
                indexList = [];
                while( b <= obj.nBlocks )
                    % get block data
                    blockX{b} = targetBlock{1,b}.data(:,1);        % same x for all targets
                    blockPDF{b} = wt(1)*targetBlock{1,b}.data(:,2);
                    blockCDF{b} = wt(1)*targetBlock{1,b}.data(:,3);

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

                %\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
                % Begin: stitch blocks together
                %/////////////////////////////////////////////////////////////////////////
                % ------------------------------ for debugging: quick look at the sections
                obj.blocks_x = blockX;
                obj.blocks_pdf = blockPDF;
                obj.blocks_cdf = blockCDF;
                obj.block_indx = indexList;

                plt_blockpdf = false;
                plt_blocksize = false;
                plt_stitchpdf = false;
                
                if plt_blockpdf

                    fig_dir = fullfile('figures_manuscript','obt_figs');

                    fig_name = 'plt_blockpdf';
                    figure('Name',fig_name)
                    hold on
                    for b=1:length(indexList)
                        plot( obj.blocks_x{obj.block_indx(b)} , obj.blocks_pdf{obj.block_indx(b)} )
                    end
                    ylabel('$\hat{f}(x)$','Interpreter','latex')
                    xlabel('$x$','Interpreter','latex')
                    if max( obj.blocks_x{obj.block_indx(length(obj.block_indx))}) < 1.1
                        ylim([0,6])
                    else
                        ylim([0,1])
                    end
                    bp = gca;
                    if save_figs
                        saveas(bp, fullfile(fig_dir, [fig_name, '.png']))
                        saveas(bp, fullfile(fig_dir, [fig_name, '.fig']))
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
                if length(indexList)-1 > 0
                    for b=1:length(indexList)-1

                        xmin = min(blockX{indexList(b+1)});
                        xmax = max(blockX{indexList(b)});
                        xmin_list = mink(blockX{indexList(b+1)},5);
                        xmax_list = maxk(blockX{indexList(b)},5);

                        %                 Lnot = or( (obj.sx <= xmin) , (obj.sx >= xmax) );
                        Lnot = or( (obj.sx < xmin) , (obj.sx > xmax) );
                        %                 Lnot = or( (obj.sx < xmin*(1+0.0001)) , (obj.sx > xmax*(1-0.0001)) );  % works for uniform, unform-mix, GP
                        %                 Lnot = or( (obj.sx < xmin*(1+0.0000001)) , (obj.sx > xmax*(1-0.0000001)) );  % test for stable
                        %                 Lnot = or( (obj.sx < xmin_list(2)) , (obj.sx > xmax_list(2)) );
                        xStitch = obj.sx(~Lnot);                          % within overlap region
                        k0 = length(xStitch);
                        if( k0 < 1 ) %51 <------------------------------ currently very small!

                            randColor = rand(length(indexList),3);
                            figure('Name','output blocks')
                            hold on;
                            for nb=1:length(indexList)
                                block_min = min(blockX{indexList(nb)});
                                block_max = max(blockX{indexList(nb)});
                                disp(['original block: ',num2str(nb), ' (',num2str(block_min),', ',num2str(block_max),')'])
                                %                         plot([block_min,block_max],[2*nb,2*nb],'-o','Color',randColor(nb,:))
                                plot([block_min,block_max],[1*nb,1*nb],'-o','Color',randColor(nb,:))

                            end


                            figure('Name','input scaled blocks')
                            hold on;
                            for nb=1:length(indexList)
                                block_dat = MatPDFsample{indexList(nb)};
                                block_min = min(block_dat);
                                block_max = max(block_dat);
                                disp(['original block: ',num2str(nb), ' (',num2str(block_min),', ',num2str(block_max),')'])
                                %                         plot([block_min,block_max],[2*nb,2*nb],'-o','Color',randColor(nb,:))
                                plot([block_min,block_max],[1*nb,1*nb],'-o','Color',randColor(nb,:))

                            end

                            figure('Name','input not scaled blocks')
                            hold on;
                            for nb=1:length(indexList)
                                block_dat = MatPDFsample_notscaled{indexList(nb)};
                                block_min = min(block_dat);
                                block_max = max(block_dat);
                                disp(['original block: ',num2str(nb), ' (',num2str(block_min),', ',num2str(block_max),')'])
                                %                         plot([block_min,block_max],[2*nb,2*nb],'-o','Color',randColor(nb,:))
                                plot([block_min,block_max],[1*nb,1*nb],'-o','Color',randColor(nb,:))

                            end

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
                        %                 u0 = min(ub);
                        %                 u0 = min(ub) + 1e-10; % ------------------------------------------------------------ TEMP FIX
                        %                 u1 = max(ub);
                        %                 v0 = min(vb);
                        %                 v0 = min(vb) + 1e-10;
                        %                 v1 = max(vb);
                        %
                        %                 u = (ub - u0)/(u1 - u0);
                        %                 v = (vb - v0)/(v1 - v0);
                        %                 power = 2;
                        %                 aLower = (1 - u).^(power);
                        %                 aUpper = v.^(power);
                        %                 temp = aLower + aUpper;
                        %
                        %                 f_lower = aLower./temp;
                        %                 f_upper = aUpper./temp;
                        %
                        %                 stitchPDF = PDFlower.*f_lower + PDFupper.*f_upper;



                        power=2;
                        u = 0: (1/(length(xStitch)-1)) : 1;
                        aLower = (1 - u).^(power);
                        aUpper = u.^(power);
                        f_lower = aLower./(aLower + aUpper);
                        f_upper = aUpper./(aLower + aUpper);
                        stitchPDF = PDFlower.*f_lower + PDFupper.*f_upper;



                        if sum(~isfinite(PDFlower)) || sum(~isfinite(obj.sPDF))|| sum(~isfinite(PDFupper)) || sum(~isfinite(stitchPDF))

                            % debugging
                            disp(['    left block # = ',num2str(b)]);
                            disp(['   right block # = ',num2str(b+1)]);
                            disp(['            xmax = ',num2str(xmax)]);
                            disp(['            xmin = ',num2str(xmin)]);
                            disp(['         overlap = ',num2str(k0)]);

                            disp(['xb: min ',num2str(min(xb)),' max ',num2str(max(xb))])
                            disp(['xStitch: min ',num2str(min(xStitch)),' max ',num2str(max(xStitch))])
                            disp(['xb1: min ',num2str(min(xb1)),' max ',num2str(max(xb1))])
                            disp(['block: ',num2str(indexList(b))])

                            ubtest = sum(~isfinite(ub))
                            vbtest = sum(~isfinite(vb))

                            ub_vals = ub(~isfinite(ub));
                            vb_vals = ub(~isfinite(vb));
                            ub_vals2 = ~isfinite(ub);
                            vb_vals2 = ~isfinite(vb);

                            randColor = rand(length(indexList),3);

                            figure('Name','input scaled blocks')
                            hold on;
                            for nb=1:length(indexList)
                                block_dat = MatPDFsample{indexList(nb)};
                                block_min = min(block_dat);
                                block_max = max(block_dat);
                                disp(['original block: ',num2str(nb), ' (',num2str(block_min),', ',num2str(block_max),')'])
                                %                         plot([block_min,block_max],[2*nb,2*nb],'-o','Color',randColor(nb,:))
                                plot([block_min,block_max],[1*nb,1*nb],'-o','Color',randColor(nb,:))

                            end

                            figure('Name','input not scaled blocks')
                            hold on;
                            for nb=1:length(indexList)
                                block_dat = MatPDFsample_notscaled{indexList(nb)};
                                block_min = min(block_dat);
                                block_max = max(block_dat);
                                disp(['original block: ',num2str(nb), ' (',num2str(block_min),', ',num2str(block_max),')'])
                                %                         plot([block_min,block_max],[2*nb,2*nb],'-o','Color',randColor(nb,:))
                                plot([block_min,block_max],[1*nb,1*nb],'-o','Color',randColor(nb,:))

                            end

                            figure('Name','output blocks')
                            hold on;
                            for nb=1:length(indexList)
                                block_min = min(blockX{indexList(nb)});
                                block_max = max(blockX{indexList(nb)});
                                %                         plot([block_min,block_max],[2*nb,2*nb],'-o','Color',randColor(nb,:))
                                plot([block_min,block_max],[1*nb,1*nb],'-o','Color',randColor(nb,:))

                                disp(['estimate block: ',num2str(nb), ' (',num2str(block_min),', ',num2str(block_max),')'])
                                x_min = min(blockX{indexList(nb+1)});
                                x_max = max(blockX{indexList(nb)});

                                %                         Lnot = or( (obj.sx <= x_min) , (obj.sx >= x_max) );
                                Lnot = or( (obj.sx < x_min) , (obj.sx > x_max) );
                                xStitch2 = obj.sx(~Lnot);
                                xStitch2_max = max(xStitch2);
                                xStitch2_min = min(xStitch2);
                                disp(['estimate stitch block: ',num2str(nb), ' (',num2str(xStitch2_min),', ',num2str(xStitch2_max),')'])
                                %                         plot([xStitch2_min,xStitch2_max],[2*(nb+1/2),2*(nb+1/2)],'-s','Color',randColor(nb,:))
                                try
                                    plot([xStitch2_min,xStitch2_max],[1*(nb+1/2),1*(nb+1/2)],'-s','Color',randColor(nb,:))
                                catch
                                    disp('temp')
                                end

                            end
                            test1 = sum(~isfinite(PDFlower))
                            test2 = sum(~isfinite(f_lower))
                            test2 = sum(~isfinite(PDFlower))
                            test3 = sum(~isfinite(f_upper))

                            test4 = sum(~isfinite(obj.sPDF))
                            test5 = sum(~isfinite(stitchPDF))
                            warning('non-finite values')

                            disp('test')
                        end


                        if plt_stitchpdf
                            fig_dir = fullfile('figures_manuscript','obt_figs');
                            publicationQuality();

                            fig_name = ['cdf_pdf_stitching_b_',num2str(b)];
                            figure('Name',fig_name)
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
                            bp = gca;
                            if save_figs
                                saveas(bp, fullfile(fig_dir, [fig_name, '.png']))
                                saveas(bp, fullfile(fig_dir, [fig_name, '.fig']))
                            end

                            set(0,'DefaultFigureColor','white')
%                             fig.InvertHardcopy = 'off';
                            width = 2;                                                                 % Width in inches
                            height = 4;                                                                % Height in inches
                            alw = 1.5;                                                                 % AxesLineWidth 
                            fsz = 14;                                                                  % Fontsize 
                            lw = 1.5;                                                                  % LineWidth 
                            msz = 8;                                                                   % MarkerSize 
                            set(0,'defaultAxesFontSize',fsz); 
                            set(0,'defaultLineLineWidth',lw);   
                            set(0,'defaultLineMarkerSize',msz); 
                            set(0,'defaultAxesLineWidth',alw);
                            defpos = get(0,'defaultFigurePosition');
                            set(0,'defaultFigurePosition', [defpos(1) defpos(2) width*100, height*100]); 
                            set(0,'defaultFigurePosition', [400, 50, width*100, height*110]); 

                            fig_name = ['cdf_pdf_stitching_cdf_b_',num2str(b)];
                            figure('Name',fig_name)
                            plot(xStitch,f_lower,'r');
                            hold on;
                            plot(xStitch,f_upper,'b');
                            ylabel('$\hat{F}(x)$','Interpreter','latex')
                            xlabel('$x$','Interpreter','latex')
                            bp = gca;
                            if save_figs
                                saveas(bp, fullfile(fig_dir, [fig_name, '.png']))
                                saveas(bp, fullfile(fig_dir, [fig_name, '.fig']))
                            end
                            %----------------------------------------------
                            fig_name = ['cdf_pdf_stitching_stitches_b_',num2str(b)];
                            figure('Name',fig_name)
                            plot(xStitch,PDFlower,'r');
                            hold on;
                            plot(xStitch,PDFupper,'b');
                            plot(xStitch,stitchPDF,'k');
                            ylabel('$\hat{f}(x)$','Interpreter','latex')
                            xlabel('$x$','Interpreter','latex')
                            bp = gca;
                            if save_figs
                                saveas(bp, fullfile(fig_dir, [fig_name, '.png']))
                                saveas(bp, fullfile(fig_dir, [fig_name, '.fig']))
                            end
                            %----------------------------------------------
                            fig_name = ['cdf_pdf_stitching_stitched_b_',num2str(b)];
                            figure('Name',fig_name)
                            plot(xStitch,stitchPDF,'k');
                            ylabel('$\hat{f}(x)$','Interpreter','latex')
                            xlabel('$x$','Interpreter','latex')
                            bp = gca;
                            if save_figs
                                saveas(bp, fullfile(fig_dir, [fig_name, '.png']))
                                saveas(bp, fullfile(fig_dir, [fig_name, '.fig']))
                            end


                            set(0,'DefaultFigureColor','white')
                            %                             fig.InvertHardcopy = 'off';
                            width = 4;                                                                 % Width in inches
                            height = 4;                                                                % Height in inches
                            alw = 1.5;                                                                 % AxesLineWidth
                            fsz = 14;                                                                  % Fontsize
                            lw = 1.5;                                                                  % LineWidth
                            msz = 8;                                                                   % MarkerSize
                            set(0,'defaultAxesFontSize',fsz);
                            set(0,'defaultLineLineWidth',lw);
                            set(0,'defaultLineMarkerSize',msz);
                            set(0,'defaultAxesLineWidth',alw);
                            defpos = get(0,'defaultFigurePosition');
                            set(0,'defaultFigurePosition', [defpos(1) defpos(2) width*100, height*100]);
                            set(0,'defaultFigurePosition', [400, 50, width*100, height*110]);


                            %----------------------------------------------
                            fig_name = ['cdf_pdf_stitching_stitches_wide_b_',num2str(b)];
                            figure('Name',fig_name)
                            plot(xStitch,PDFlower,'r');
                            hold on;
                            plot(xStitch,PDFupper,'b');
                            plot(xStitch,stitchPDF,'k');
                            ylabel('$\hat{f}(x)$','Interpreter','latex')
                            xlabel('$x$','Interpreter','latex')
                            bp = gca;
                            if save_figs
                                saveas(bp, fullfile(fig_dir, [fig_name, '.png']))
                                saveas(bp, fullfile(fig_dir, [fig_name, '.fig']))
                            end


                            publicationQuality();


                        end


                        [~,indx] = min( abs(obj.sx-xmin) );
                        dk = k0 - 1;
                        obj.sPDF(indx:indx+dk) = stitchPDF;

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

                    end
                else
%                         [~,indx] = min( abs(obj.sx-xmin) );
%                             dk = k0 - 1;
%                             obj.sPDF(indx:indx+dk) = stitchPDF;

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
                        obj.sCDF = pNorm*(obj.sCDF/temp) + pL;
                end
            else

                try
                    [~, obj.sx, obj.sPDF, obj.sCDF, ~,lagrange] = EstimatePDF(inputSample);
                    obj.failed = 0;
                catch
                    warning('Problem using function.  Assigning a value of 0.');
                    lagrange = 0;
                    obj.sx = linspace(min(inputSample),max(inputSample),100);
                    obj.sPDF = 0*ones(100,1);
                    obj.sCDF = 0*ones(100,1);
                    obj.failed = 1;
                end
                %                 LG(1,b) = size(lagrange,1);
                %                 LG_vals{b} = lagrange;

                obj.LG_vals{1} = lagrange;
                obj.LG = [size(lagrange,1)];
                obj.LG_max = max(obj.LG);
                obj.LG_sum = sum(obj.LG);

                obj.nBlocks = 1;
                obj.blocks_x{1} = obj.sx;
                obj.blocks_pdf{1} = obj.sPDF;
                obj.blocks_cdf{1} = obj.sCDF;
                obj.block_indx = [1];
                obj.block_size = obj.N;
                obj.block_scale = 0.5*( obj.N - 1 );

            end



            % calculate information to plot      section 19
            sample = x(2:end-1);
            % calculate & plot SQR

            % adjust u range
            sampleUpLim = max(obj.sx);
            sampleLoLim = min(obj.sx);
            [row, ~] = find(sample <= sampleUpLim & sample >= sampleLoLim);
            try
                obj.u = interp1(obj.sx,obj.sCDF,sample(row));    % get corresponding u for each x in sample
            catch
                test ='test';
            end
            uref = (1:size(sample(row),1))/(size(sample(row),1) - 1);
            if( size(uref,1) ~= size(obj.u,1) )
                obj.u = obj.u';
            end


            % get scaled residual
            obj.sqr = sqrt(obj.N)*(obj.u - uref); % normal formula has sqrt(N+2) but N -> N-2
            [obj.u,obj.sqr] = utils.sqr(obj.sx,obj.sPDF,inputSample);
        end
    end
end