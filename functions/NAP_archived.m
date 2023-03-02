classdef NAP
    properties (Constant)
        % hyper parameter set
        p1 = 1;
        p2 = 0.55;
%         p2 = 1/2;
        p3 = 1;
        p4 = 0.33;
%         p4 = 1/3;
        p5 = 2;
%         p6 = 0.0625;
        p6 = 1/8;
        p7 = 0.5;
        p8 = 40;
    end
    properties (SetAccess = public)
        % maximum allowed block size
        max_bs = 1e5;
    end
    properties
        % estimate data
        failed;
        sx;
        pdf;
        cdf;
        u;
        sqr;
        % possible outputs
        p = [NAP.p1,NAP.p2,NAP.p3,NAP.p4,NAP.p5,NAP.p6,NAP.p7,NAP.p8];
        n_blocks;
        N;
        binN;
        T;
        xi_lvl;
        xi0;
        % estimate data per block
        blocks_x;
        blocks_pdf;
        blocks_cdf;
        block_indx;
        % meta data
        lagrange_vals;
        lagrange;
        lagrange_max;
        lagrange_sum;
        block_size;
        block_scale;
    end
    methods
        function [obj] = stitch(obj, sample_origin, serial)

            % initialze the obj.failed trip flag to false
            obj.failed = 0;

            % user input
            trgt_coverage = round(70);

            % minimum number of allowed sample size
            ns_min = 50; 

            % some error checking
            trgt_coverage = sort(trgt_coverage);
            n_trgts = length(trgt_coverage);

            % universal scoring function
            loglikeli = utils.likelihood();
            x0_score = loglikeli(:,1);
            y0_coverage = loglikeli(:,2);
            x_score = x0_score;
            y_coverage = 100*y0_coverage;

            % define weights
            trgt_bound_low = trgt_coverage - 2.5;
            trgt_bound_high = trgt_coverage + 2.5;
            if( trgt_bound_low < 0.099999 )
                error('A target coverage value is too low: Lowest = 2.6');
            end
            if( trgt_bound_high > 95.000001 )
                error('A target coverage value is too high: Highest = 92.5');
            end

            % assign weights
            idx = 15:length(x_score)-15;
            deltax = x_score(idx+14) - x_score(idx-14);
            deltay = y_coverage(idx+14) - y_coverage(idx-14);
            wScore = deltay./deltax;
            wt = size(trgt_coverage,1);
            for k=1:n_trgts
                a = y_coverage - trgt_coverage(k);
                [~,idx] = min( abs(a) );
                wt(k) = wScore(idx);
            end
            wtsum = sum(wt);
            wt = wt/wtsum;

            % check sample for unique values
            sample = unique(sample_origin);
            n_duplicates = length(sample_origin) - length(sample);
            if n_duplicates > 0
                warning(['Number of dupliates: ', num2str(2*n_duplicates)])
            end

            % sort sample and delta-x array
            s_sorted = sort(sample);
            dx = sample(2:end) - sample(1:end-1);
            dxs = sort(dx);

            obj.N = length(s_sorted);
            if( obj.N < ns_min )
                error(['Sample size is too small. ns_min = ',num2str(ns_min)]);
            end
            obj.binN = 15;
%             obj.binN = 300;

            % OBSERVATION WINDOW: CDF Boundary Probabilities
            % probability for data to be  left of window
            prob_left = 0.5/obj.N;
            % probability for data to be right of window
            prob_right = 0.5/obj.N;
            % probability for data to fall within window
            prob_norm = 1 - prob_left - prob_right;

            % partition data into primary blocks
            obt_blocks = blocks;
            obt_blocks.sample = s_sorted;
            obt_blocks.dx = dx;
            obt_blocks.dxs = dxs;
            obt_blocks.binNs = obj.binN;
            [j,obj.n_blocks,p_list,obj.T, obj.xi_lvl,obj.xi0] ...
                = obt_blocks.bin_width_size(obt_blocks);

            partition_low = zeros(1,obj.n_blocks);
            partition_high = zeros(1,obj.n_blocks);
            % create staggered (secondary) blocks
            if( obj.n_blocks > 2 )
                % block size of three is a special case
                if( obj.n_blocks == 3 )
                    partition_low(1) = 1;
                    partition_high(1) = p_list(1);
                    partition_low(2) = round(0.5*(1  + p_list(1)));
                    partition_high(2) = round(0.5*(p_list(1) + obj.N));
                    partition_low(3) = p_list(1);
                    partition_high(3) = obj.N;
                else
                    % set index for staggered blocks (sb)
                    sb_idx = 1;
                    % set edge positions for first sb
                    partition_low(sb_idx) = 1;
                    partition_high(sb_idx) = p_list(1);
                    sb_idx = sb_idx + 1;
                    % set edge positions for second sb as average first 
                    % blocks (fb)
                    partition_low(sb_idx) = ...
                        round(0.5*(1 + p_list(1)));
                    partition_high(sb_idx) = ...
                        round(0.5*(p_list(1) + p_list(2)));

                    % Loop to decide other sb edge positions
                    for i=1:j-2
                        % update indexer for kBlock indices
                        sb_idx = sb_idx + 1;
                        % set edge positions as first block for i-th sb
                        partition_low(sb_idx) = p_list(i);
                        partition_high(sb_idx) = p_list(i+1);
                        % update indexer for kBlock indices
                        sb_idx = sb_idx + 1;
                        % set edge positions as average of fb positions
                        % for i+1-th sb
                        partition_low(sb_idx) = ...
                            round(0.5*(p_list(i) + p_list(i+1)));
                        partition_high(sb_idx) = ...
                            round(0.5*( p_list(i+1) + p_list(i+2)));
                    end
                    sb_idx = sb_idx + 1;
                    partition_low(sb_idx) = p_list(j-1);
                    partition_high(sb_idx) = p_list(j);
                    sb_idx = sb_idx + 1;
                    partition_low(sb_idx) = ...
                        round(0.5*(p_list(j-1) + p_list(j)));
                    partition_high(sb_idx) = ...
                        round(0.5*(p_list(j) + obj.N ));

                    partition_low(obj.n_blocks) = partition_high(end-2);
                    partition_high(obj.n_blocks) = obj.N;
                end
            else % => block size is 1
                partition_low(1) = 1;
                partition_high(1) = obj.N;
            end

            % TODO: update
            % Determine block sizes
            block_size_list = partition_high - partition_low + 1;
            flag = 1;
            if( block_size_list(1) < obj.binN )
                flag = -1;
                k = obj.binN - block_size_list(1);
                tt = block_size_list - k;
                for idx1=2:obj.n_blocks
                    if( tt(idx1) > obj.binN )
                        idx_ref = idx1;
                        if( idx_ref == obj.n_blocks )
                            error(['All blocks are too small! ',...
                                'Sample is too hard']);
                        end
                        break;
                    end
                end
                partition_low(1) = 1;
                partition_high(1) = partition_high(1) + k;
                for idx2=2:idx_ref-1
                    partition_high(idx2) = partition_high(idx2) + k;
                    partition_low(idx2) = partition_low(idx2) + k;
                end
                partition_low(idx_ref) = partition_low(idx_ref) + k;
                partition_low(idx_ref+1) = partition_low(idx_ref+1) + k;
            end
            if( block_size_list(obj.n_blocks) < obj.binN )
                flag = -1;
                k = obj.binN - block_size_list(obj.n_blocks);
                tt = block_size_list - k;
                for idx3=obj.n_blocks-1:-1:1
                    if( tt(idx3) > obj.binN )
                        idx_ref = idx3;
                        if( idx_ref == 1 )
                            error(['All blocks are too small! ',...
                                'Sample is too hard']);
                        end
                        break;
                    end
                end
                partition_high(obj.n_blocks) = obj.N;
                partition_low(obj.n_blocks) =...
                    partition_low(obj.n_blocks) - k;
                for idx4=obj.n_blocks-1:-1:idx_ref+1
                    partition_high(idx4) = partition_high(idx4) - k;
                    partition_low(idx4) = partition_low(idx4) - k;
                end
                partition_high(idx_ref) = partition_high(idx_ref) - k;
                partition_high(idx_ref-1) = partition_high(idx_ref-1) - k;
            end
            if( flag < 0 )
                block_size_list = partition_high - partition_low + 1;
            end

            % get length scale & shifts for each block

            % length of s_sorted as the span of a block
            blockScale = zeros(1,obj.n_blocks);
            % shift in s_sorted to center the span about 0
            blockShift = zeros(1,obj.n_blocks);
            for b=1:obj.n_blocks
                p_low = partition_low(b);
                p_high = partition_high(b);
                blockScale(b) = 0.5*( s_sorted(p_high) - s_sorted(p_low) );
                blockShift(b) = 0.5*( s_sorted(p_high) + s_sorted(p_low) );
            end

            % track meta data
            obj.block_size = block_size_list;
            obj.block_scale = blockScale;

            % set up probability factors
            nBs0 = block_size_list - 1;
            nBs0(1) = nBs0(1) + 0.5;
            nBs0(obj.n_blocks) = nBs0(obj.n_blocks) + 0.5;

            % initialize tempStruc for parfor and generate block
            % sub-samples
            bounds = cell(1,obj.n_blocks);
            subsamples = cell(1,obj.n_blocks);
            for b_idx1=1:obj.n_blocks
                j1 = partition_low(b_idx1);
                j2 = partition_high(b_idx1);
                s = 1/blockScale(b_idx1);

                subsamples{b_idx1} = ...
                    s*( s_sorted(j1+1:j2-1) - blockShift(b_idx1) )';
                
                % bound block estiamtes for interior/exterior blocks
                if b_idx1 == 1
                    tempStruc.highBound = 1;
                elseif b_idx1 == obj.n_blocks
                    tempStruc.lowBound = -1;
                else
                    tempStruc.lowBound = -1;
                    tempStruc.highBound = 1;
                end
                bounds{b_idx1} = tempStruc;
            end
            
            % initialize vector to hold all lagrainge mutiplers per block
            lagrange = zeros(1,obj.n_blocks);
            lagrange_vals = cell(1,obj.n_blocks);

            % conditionally execute parfor or for loop
            if serial
                parforArg = 0;
            else
                parforArg = Inf;
            end
            parfor (b=1:obj.n_blocks, parforArg)
                for t=1:n_trgts
                    try
                        [~, targetBlock{t,b}.data(:,1), ...
                            targetBlock{t,b}.data(:,2), ...
                            targetBlock{t,b}.data(:,3), ~,lagrange] = ...
                            EstimatePDF(subsamples{b}, bounds{b});    
                    catch
                        warning(['Problem using function.',...
                            ' Assigning a value of 0.',' t: ',...
                            num2str(t),' b: ',num2str(b)]);
                        lagrange = 0;
                        targetBlock{t,b}.data(:,1) = ...
                            linspace(min(subsamples{b}),...
                                max(subsamples{b}),100);
                        targetBlock{t,b}.data(:,2) = 0*ones(100,1);
                        targetBlock{t,b}.data(:,3) = 0*ones(100,1);
                    end
                    lagrange(1,b) = size(lagrange,1);
                    lagrange_vals{b} = lagrange;
                end
            end

            % track meta data
            obj.lagrange_vals = lagrange_vals;
            obj.lagrange = lagrange;
            obj.lagrange_max = max(lagrange);
            obj.lagrange_sum = sum(lagrange);

            % combine different estimate data for all blocks
            b_pdf = cell(1,obj.n_blocks);
            b_cdf = cell(1,obj.n_blocks);
            b_x = cell(1,obj.n_blocks);

            b = 1;
            b_idx2 = [];
            while( b <= obj.n_blocks )
                % get block data
                b_x{b} = targetBlock{1,b}.data(:,1);
                b_pdf{b} = wt(1)*targetBlock{1,b}.data(:,2);
                b_cdf{b} = wt(1)*targetBlock{1,b}.data(:,3);

                for t=2:n_trgts
                    b_pdf{b} = b_pdf{b} + wt(t)*targetBlock{t,b}.data(:,2);
                    b_cdf{b} = b_cdf{b} + wt(t)*targetBlock{t,b}.data(:,3);
                end

                % scale block data

                % Note:  sInv = 1/s = 1/blockScale(b)
                sInv = blockScale(b);
                s = 1/sInv;
                % Wb = nBs0/N = weight based on frequency count
                sWb = s*(nBs0(b)/obj.N);

                % apply scaling
                b_pdf{b} = sWb*b_pdf{b};
                b_x{b} = sInv*b_x{b} + blockShift(b);

                % b_idx2 is used to evaluate blocks later on
                b_idx2 = [b_idx2,b];
                b = b + 1;
            end

            obj.n_blocks = b-1;

            %\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
            % Begin: stitch blocks together

            obj.blocks_x = b_x;
            obj.blocks_pdf = b_pdf;
            obj.blocks_cdf = b_cdf;
            obj.block_indx = b_idx2;

            % sync block data into a single range
            obj.sx = [];
            obj.pdf = [];
            for b=1:length(b_idx2)
                obj.sx = [obj.sx, b_x{b_idx2(b)}'];
                obj.pdf = [obj.pdf, b_pdf{b_idx2(b)}'];
            end
            [obj.sx,idx] = unique(obj.sx);
            obj.pdf = obj.pdf(idx);

            % stitch the blocks together
            for b=1:length(b_idx2)-1

                xmin = min(b_x{b_idx2(b+1)});
                xmax = max(b_x{b_idx2(b)});
                Lnot = or( (obj.sx < xmin) , (obj.sx > xmax) );
                % within overlap region
                x_stitch = obj.sx(~Lnot); 
                k0 = length(x_stitch);

                % overlap is too small or does not overlap
                if( k0 < 10 )
                    disp(['    left block # = ',num2str(b)]);
                    disp(['   right block # = ',num2str(b+1)]);
                    disp(['            xmax = ',num2str(xmax)]);
                    disp(['            xmin = ',num2str(xmin)]);
                    disp(['         overlap = ',num2str(k0)]);
                    warning('overlap is too small!');
                    obj.failed = 1;
                end
                
                % interpolate PDFs over stitched region
                [xb,idx] = unique(b_x{b_idx2(b)});
                yb = b_pdf{b_idx2(b)}(idx);
                PDFlower = interp1(xb,yb,x_stitch);

                [xb1,idx] = unique(b_x{b_idx2(b+1)});
                yb1 = b_pdf{b_idx2(b+1)}(idx);
                PDFupper = interp1(xb1,yb1,x_stitch);

                % stitch PDFs
                power=2;
                u = 0: (1/(length(x_stitch)-1)) : 1;
                aLower = (1 - u).^(power);
                aUpper = u.^(power);
                f_lower = aLower./aLower + aUpper;
                f_upper = aUpper./aLower + aUpper;
                pdf_stitch = PDFlower.*f_lower + PDFupper.*f_upper;

                [~,idx] = min( abs(obj.sx-xmin) );
                dk = k0 - 1;
                obj.pdf(idx:idx+dk) = pdf_stitch;

            end

            % calcuate CDF
            obj.cdf = zeros( size(obj.pdf) );
            obj.cdf(1) = 0;
            kmax = length(obj.cdf);
            for k=2:kmax
                fave = 0.5*( obj.pdf(k) + obj.pdf(k-1) );
                area = fave*( obj.sx(k) - obj.sx(k-1) );
                obj.cdf(k) = obj.cdf(k-1) + area;
            end
            % recalling what prob_left, prob_right and prob_norm are
            obj.cdf = prob_norm*(obj.cdf/obj.cdf(kmax)) + prob_left;

            % END: stitch blocks together
            %//////////////////////////////////////////////////////////////

            % calculate SQR
            s_sqr = s_sorted(2:end-1);
            [row, ~] = find(s_sqr <= max(obj.sx) & s_sqr >= min(obj.sx));
            % get corresponding u for each s_sorted in sample
            obj.u = interp1(obj.sx,obj.cdf,s_sqr(row));
            uref = (1:size(s_sqr(row),1))/(size(s_sqr(row),1) - 1);
            if( size(uref,1) ~= size(obj.u,1) )
                obj.u = obj.u';
            end

            % get scaled residual
            % normal formula has sqrt(N+2) but N -> N-2
            obj.sqr = sqrt(obj.N)*(obj.u - uref);
            [obj.u,obj.sqr] = utils.sqr(obj.sx,obj.pdf,s_sqr);
        end
    end
end