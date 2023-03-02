classdef blocks < NAP % inherit properties i.e. p vector and max block size

    properties
        sample;
        dx;
        dxs;
        Ns;
        window;
        min_bin;
        n_bin;
    end
    methods(Static)
        function [j,n_blocks,p_list,gamma,xi_lvl,xi0] = ...
                bin_width_size(obj)

            obj.Ns = length(obj.sample);

            if( obj.Ns < obj.n_bin )
                n_blocks = 1;
            else
                % determine number of blocks and their sizes
                [p_list,gamma,xi_lvl,xi0] = obj.r_tree(obj);
                % remove tiny left and right
                if( p_list(1) < obj.n_bin )
                    p_list = p_list(2:end);
                end
                if( obj.Ns - p_list(end) + 1 < obj.n_bin )
                    p_list = p_list(1:end-1);
                end
                j = length(p_list);
                n_blocks = 2*j + 1;
            end
        end

        function [plist,gamma,xi_lvl,xi0] = r_tree(obj)

            % create gamma threshold
            gamma = obj.p(5)*obj.Ns^obj.p(4);
            obj.window = ceil(obj.p(6)*obj.Ns^obj.p(7));
            maxLevel = obj.p(8);

            % minimum blocksize
            obj.min_bin = ceil(2*obj.window);
            if obj.min_bin > obj.Ns
                error('sample size too small for window size')
            end

            % track number of branches per level
            nbranch = 1;

            % initialize vector to track all created partitions

            % set end points of sample length as partition left (pL) and
            % partion right (pR)
            plist = [1 obj.Ns];
            % track every attempted partition for all levels
            plevel = {{[1; obj.Ns]}};
            % initialize array to track newly created partitions
            % for plotting purposes
            pdiff = {[1; obj.Ns]};
            % calcualte inital BR of intire sample
            B0 = obj.Ns;
            R0 = obj.get_ratio(obj, obj.sample);
            xi0 = obj.br_product(obj, B0, R0, obj.Ns);

            % clear array to hold all BR values per block per level
            xi_lvl = {xi0};

            % beggin level loop
            if xi0 > gamma
                for jj = 1:maxLevel
                    % vector to hold all attempted partitions per level
                    p_lvl_hold = [];
                    % vector to hold all BR values per level
                    xi_hold = [];
                    % beggin branch loop
                    for b = 1:nbranch
                        % define block size (B)
                        B = plist(b+1) - plist(b);
                        if obj.min_bin >= B
                            continue
                        end
                        % update left
                        boundry_left = obj.min_bin + 1;
                        boundry_right = B - obj.min_bin;

                        % block to small for for minimization
                        % given window size
                        if boundry_right - boundry_left < 3
                            break;
                        end

                        % golden ration bifraction minimization
                        [xi_left,xi_right,partition] =...
                            obj.min_br_gold(obj,...
                            obj.sample(plist(b):plist(b+1)));

                        % find minimum BR for two newly created blocks
                        BR = min(xi_left,xi_right);

                        % update block boundaries of sample with correctly
                        % placed partition
                        newPar = plist(b) + partition;

                        % hold all BRs per level for later evaluation
                        xi_hold = [xi_hold,BR];

                        % un-balanced tree
                        if BR >= gamma
                            p_lvl_hold = [p_lvl_hold, newPar];
                        end
                    end
                    
                    plist = sort([plist p_lvl_hold]);
                    % update nbranch
                    nbranch = length(plist) - 1;

                    % assign partition list to array for plotting
                    plevel{jj+1,1} = {plist'};

                    % exit for special cases where B < min_bin
                    % or b_right - b_left < 3
                    if isempty(xi_hold)
                        break;
                    end
                    % assign BR per level to array for plotting
                    xi_lvl{jj+1,1} = xi_hold;

                    % find newly accepted partitions
                    [C,~] = setdiff(plevel{jj+1,1}{1,1}(:,1),...
                        plevel{jj,1}{1,1}(:,1));

                    % STOP SEARCH: if no new partitions are accepted
                    if isempty(C)
                        break;
                    end

                    % update changes with newly created partitions
                    pdiff{jj+1,1} = C;
                end
                plist = plist';

                % SPLITTING ROUTINE FOR LARGE SAMPLES

                % vector to to add new partitons too
                lrg_n_check = plist;

                % vector to hold updated partition list
                holder = plist;

                % while loop flag
                run_split = true;

                while run_split
                    % triggers exit flag for while loop
                    split_count = 0;

                    % loop over modified partition list (holder)
                    for k = 1:length(holder)-1

                        % add partion between elements when diff > max_bs
                        if holder(k+1)-holder(k) > obj.max_bs

                            split = floor((holder(k+1)-holder(k))/2);
                            % update new partiton list

                            lrg_n_check = [lrg_n_check;...
                                lrg_n_check(k)+ split];
                            % update counter: number of found splits

                            split_count = split_count + 1;
                        end
                    end
                    lrg_n_check = sort(lrg_n_check);
                    holder = lrg_n_check;

                    % if no splits exit routine
                    if split_count == 0
                        run_split = false;
                    end
                end
                plist = lrg_n_check;
            end
        end

        function [xi_left,xi_right,partition] = min_br_gold(obj, sample)

            Ns = length(sample);
            b_left = 1 + obj.min_bin;
            b_right = Ns - obj.min_bin;
            gold_ratio = (1+sqrt(5))/2;
            dxbr = [];
            loop_count = 0;
            partition = ceil(Ns/2);

            while b_right-b_left > 2

                % define partitions
                b_left1 = partition;
                b_right1 = Ns - partition;
                r_right1 = obj.get_ratio(obj, sample(partition:end));
                r_left1 = obj.get_ratio(obj, sample(1:partition));
                dxbr_center = abs(obj.br_product(obj, b_left1, r_left1, Ns)...
                    -obj.br_product(obj, b_right1, r_right1, Ns));
                % --------------------
                p_left = partition - 1;
                b_left2 = p_left;
                b_right2 = Ns - p_left ;
                r_right2 = obj.get_ratio(obj, sample(p_left:end));
                r_left2 = obj.get_ratio(obj, sample(1:p_left));
                dxbr_left = abs(obj.br_product(obj, b_left2,r_left2,Ns)...
                    -obj.br_product(obj, b_right2,r_right2,Ns));
                % --------------------
                p_right = partition + 1;
                b_left3 = p_right;
                b_right3 = Ns - p_right;
                r_right3 = obj.get_ratio(obj, sample(p_right:end));
                r_left3 = obj.get_ratio(obj, sample(1:p_right));
                dxbr_right = abs(obj.br_product(obj, b_left3,r_left3,Ns)...
                    -obj.br_product(obj, b_right3,r_right3,Ns));
                % --------------------
                dx_cr = dxbr_center - dxbr_right;
                dx_lc = dxbr_left - dxbr_center;
                dx_lr = dxbr_left - dxbr_right;

                % parition update conditions
                % dxbr_right is smallest
                if dx_cr >= 0 && dx_lr > 0
                    b_left = p_left;
                    if b_right-b_left < 2
                        dxbr = horzcat(dxbr,[dxbr_center;b_left1;r_left1;...
                            b_right1;r_right1]);
                    else
                        % Shrink --->
                        dxbr = horzcat(dxbr,[dxbr_right;b_left3;r_left3;...
                            b_right3;r_right3]);
                        partition = round((b_left+b_right*gold_ratio)/...
                            (1+gold_ratio));
                    end
                end
                % dxbr_left is smallest
                if dx_lc <= 0 && dx_lr < 0
                    b_right = p_right;
                    if b_right-b_left < 2
                        dxbr = horzcat(dxbr,[dxbr_center;b_left1;r_left1;...
                            b_right1;r_right1]);
                    else
                        % <--- Shrink
                        dxbr = horzcat(dxbr,[dxbr_left;b_left2;r_left2;...
                            b_right2;r_right2]);
                        partition = round((b_left*gold_ratio+b_right)/...
                            (1+gold_ratio));
                    end
                end
                % dxbr_center is smallest
                if dx_lc > 0 && dx_cr < 0
                    dxbr = horzcat(dxbr,[dxbr_center;b_left1;r_left1;...
                        b_right1;r_right1]);
                    break
                end
                % dxbr_center = dxbr_left = dxbr_right = 0
                if dx_lc == 0 || dx_cr == 0 || dx_lr == 0
                    dxbr = horzcat(dxbr,[dxbr_center;b_left1;r_left1;...
                        b_right1;r_right1]);
                    break
                end
                loop_count = loop_count + 1;
            end

            dxbr = dxbr';
            xi_left = obj.br_product(obj, dxbr(end,2),dxbr(end,3),Ns);
            xi_right = obj.br_product(obj, dxbr(end,4),dxbr(end,5),Ns);
            partition = dxbr(end,2);
        end

        function r = get_ratio(obj, sample)
            % create boolean mask for subsample from entire sample
            mask = (min(sample)<=obj.sample & max(sample)>=obj.sample);

            % change one binary 1 value to 0 to create mask for difference
            % array indexes where mask is true i.e. 1
            mask_idx= find(mask==1);

            % change last true value on the right to false
            mask(mask_idx(end)) = 0;
            mask = mask(1:end-1);

            % get the dx values for this subsample in sorted order
            sort_mask = ismember(obj.dxs, obj.dx(mask)');
            dx = obj.dxs(sort_mask);

            dxMin = mean(dx(1:obj.window));
            dxMax = mean(dx(end-obj.window+1:end));
            % -------------------------------------------------------------
            % r -> inf, dxMin << dxMax and r -> 1, dxMin ~ dxMax
            % -------------------------------------------------------------
            r = dxMax/dxMin;
        end

        function BR = br_product(obj, b, r, n)
            BR = (b^obj.p(1)*r^obj.p(2))/n^obj.p(3);
        end

    end
end

