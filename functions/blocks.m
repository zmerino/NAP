classdef blocks < NAP % inherit NAP properties i.e. p vector and max block size
    properties (Constant)
        % plot anything
        plt_figs = false;
        % plot xi values per parition or number of branches per level
        plt_xi = false;
        plt_tree = false;
        save_figs = false;
    end
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
                        [brL,brR,partition] = obj.min_br_gold(obj, obj.sample(pList(b):pList(b+1)));
%                         [brL,brR,partition] = obj.gss(obj, obj.sample(pList(b):pList(b+1)));

                        % plot br and r for various paritions of whole
                        % sample
                        if jj == 1 && obj.plt_figs && obj.plt_xi
                            r_left = [];
                            r_right = [];
                            b_left = [];
                            b_right = [];
                            br_left = [];
                            br_right = [];


                            hw = waitbar(0,'Calculate Left Partions: ');

                            x_partition_left = [];
                            for n = 2*obj.window:obj.Ns-2*obj.window

                                waitbar(n /obj.Ns,hw, sprintf('Calculate Left Partions: %i/%i',n, obj.Ns));
                                B0 = n;
                                b_left = [b_left, B0];
                                R0 = obj.get_ratio(obj, obj.sample(1:B0));
                                r_left = [r_left, R0];
                                BR0 = obj.br_product(obj, B0, R0, obj.Ns);
                                br_left = [br_left, BR0];
                                x_partition_left = [x_partition_left, n];
                            end

                            hw = waitbar(0,'Calculate Right Partions: ');

                            x_partition_right = [];
                            for n = 2*obj.window:obj.Ns-2*obj.window

                                waitbar(n /obj.Ns,hw, sprintf('Calculate Right Partions: %i/%i',n, obj.Ns));

                                B0 = n;
                                b_right = [b_right, B0];
                                R0 = obj.get_ratio(obj, obj.sample(end-B0:end));
                                r_right = [r_right, R0];
                                BR0 = obj.br_product(obj, B0, R0, obj.Ns);
                                br_right = [br_right, BR0];
                                x_partition_right = [x_partition_right, obj.Ns-n];
                            end

                            figure('Name', 'R verse partition')
                            hold on;
                            plot(x_partition_left, r_left, '-b')
                            plot(x_partition_right, r_right, '-r')
                            xlabel('Partition', Interpreter='latex')
                            ylabel('R', Interpreter='latex')

                            figure('Name', 'B verse partition')
                            hold on;
                            plot(x_partition_left, b_left, '-b')
                            plot(x_partition_right, b_right, '-r')
                            xlabel('Partition', Interpreter='latex')
                            ylabel('Block Size', Interpreter='latex')

            
                            fig_dir = fullfile('figures_manuscript','obt_figs');
            
                            fig_name = 'br_verse_partition';
                            figure('Name',fig_name)
                            hold on;
                            plot(x_partition_left, br_left, '--b')
                            plot(x_partition_right, br_right, '--r')
                            plot(x_partition_right, abs(br_right-flip(br_left)), '-k')
                            xlabel('Partition', Interpreter='latex')
                            ylabel('$| \Delta \xi |$', Interpreter='latex')
                            xlim([min(x_partition_right),max(x_partition_right)])
                            bp = gca;
                            bp.YAxis.Exponent = 3;
                            bp.XAxis.Exponent = 3;
                            if obj.save_figs
                                saveas(bp, fullfile(fig_dir, [fig_name, '.png']))
                                saveas(bp, fullfile(fig_dir, [fig_name, '.fig']))
                            end
                        end

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
                pList = LargNcheck;
            end

            % FIGURES -----------------------------------------------------
            if obj.plt_figs && obj.plt_tree
                
                fig_dir = fullfile('figures_manuscript','obt_figs');


                width = 6; % Width in inches
                height = 4; 
                res = 1000;
                defpos = get(groot,'defaultFigurePosition');
                
                fig_name = 'br_values_per_level';
                f = figure('Name',fig_name);
                f.Position = [defpos(1), defpos(2), width*100, height*100];
                hold on
                plot(0:size(plevel,1)-1,log(T*ones(size(plevel,1))), '--r');
                for k = 1:size(BRlevel,1)
                    plot((k-1)*ones(size(BRlevel{k,1}(1,:),1),1), log(BRlevel{k,1}(1,:)),...
                        'o',...
                        'MarkerEdgeColor',[0,0,0],...
                        'MarkerFaceColor',[0,0,0],...
                        'MarkerSize',4,'DisplayName','none')
                    levelTrack = 1:size(plevel,1);
                end
                str = cell(1,size(levelTrack,2));
                for ii = 1:length(levelTrack)-1
                    str{ii} = sprintf('%1.0f',levelTrack(ii));
                end
                xticks(levelTrack)
                xticklabels(str)
                ylabel('ln($\xi$)', 'Interpreter','latex')
                xlabel('Tree Level', 'Interpreter','latex')
                legend('$\Gamma$', 'Interpreter','latex')
                bp = gca;
                if obj.save_figs
%                     saveas(bp, fullfile(fig_dir, [fig_name, '.png']))
                    saveas(bp, fullfile(fig_dir, [fig_name, '.fig']))
                    exportgraphics(bp,fullfile(fig_dir, [fig_name, '.png']),'Resolution',res)
                end

                % set common x limits for subplots
                sample_boundry = obj.sample(pdiff{1,1}(:,1));
                x_min = min(sample_boundry);
                x_max = max(sample_boundry);

                fig_name = 'tree_branching';
                figure('Name',fig_name)
                subplot(2,1,1)
                histogram(obj.sample)
                bp = gca;
%                 bp.YAxis.Scale ="log";
                bp.YAxis.Exponent = 3;
                ylabel('$N_b$','Interpreter','latex')
                xlim([x_min x_max])


                subplot(2,1,2)
                hold on
                % branching level track markers
                for k = 1:size(plevel,1)-1
                    plot(obj.sample(plevel{k,1}{1,1}(:,1)),...
                        (size(plevel,1)-k)*...
                        ones(size(plevel{k,1}{1,1}(:,1),1),1),...
                        'o',...
                        'MarkerEdgeColor',[0.6,0.6,0.6],...
                        'MarkerFaceColor',[0.6,0.6,0.6],...
                        'MarkerSize',5)

                    levelTrack = 0:size(plevel,1)-1;
                end
                plot(obj.sample(pList),...
                    zeros(length(pList),1),...
                    'o',...
                    'MarkerEdgeColor',[1,0,0],...
                    'MarkerFaceColor',[1,0,0],...
                    'MarkerSize',5)

                % boundries of sample markers
                plot(obj.sample(pdiff{1,1}(:,1)),...
                    (size(plevel,1)-1)*...
                    ones(size(pdiff{1,1}(:,1),1),1),...
                    'o',...
                    'MarkerEdgeColor',[0,0,0],...
                    'MarkerFaceColor',[0.6,0.6,0.6],...
                    'MarkerSize',8)
                % new partion markers
                for k = 2:size(pdiff,1)
                    plot(obj.sample(pdiff{k,1}(:,1)),...
                        (size(plevel,1)-k)*...
                        ones(size(pdiff{k,1}(:,1),1),1),...
                        'o',...
                        'MarkerEdgeColor',[0,0,0],...
                        'MarkerFaceColor',[0,0,0],...
                        'MarkerSize',8)
                end

                str = cell(1,size(levelTrack,2));
                for ii = 1:length(levelTrack)
                    str{ii} = sprintf('%1.0f',levelTrack(end+1-ii));
                end

                yticks(levelTrack)
                yticklabels(str)
                xlabel('x Range','Interpreter','latex')
                ylabel('Tree Level','Interpreter','latex')
                xlim([x_min x_max])

                bp = gca;
                if obj.save_figs
%                     saveas(bp, fullfile(fig_dir, [fig_name, '.png']))
                    saveas(bp, fullfile(fig_dir, [fig_name, '.fig']))
                    exportgraphics(bp,fullfile(fig_dir, [fig_name, '.png']),'Resolution',res)
                end


                % split figure for tree branching
                width = 6; % Width in inches
                height = 2; 
                defpos = get(groot,'defaultFigurePosition');
                
                fig_name = 'tree_branching_top';
                f = figure('Name',fig_name);
                f.Position = [defpos(1), defpos(2), width*100, height*100];
                histogram(obj.sample)
                bp = gca;
%                 bp.YAxis.Scale ="log";
                bp.YAxis.Exponent = 3;
                ylabel('$N_b$','Interpreter','latex')
                xlim([x_min x_max])

                if obj.save_figs
%                     saveas(bp, fullfile(fig_dir, [fig_name, '.png']))
                    saveas(bp, fullfile(fig_dir, [fig_name, '.fig']))
                    exportgraphics(bp,fullfile(fig_dir, [fig_name, '.png']),'Resolution',res)
                end

                fig_name = 'tree_branching_bottom';
                f = figure('Name',fig_name);
                f.Position = [defpos(1), defpos(2), width*100, height*100];
%                 set(f,'defaultFigurePosition', [defpos(1) defpos(2) width*100, height*100]); 
                hold on
                % branching level track markers
                for k = 1:size(plevel,1)-1
                    plot(obj.sample(plevel{k,1}{1,1}(:,1)),...
                        (size(plevel,1)-k)*...
                        ones(size(plevel{k,1}{1,1}(:,1),1),1),...
                        'o',...
                        'MarkerEdgeColor',[0.6,0.6,0.6],...
                        'MarkerFaceColor',[0.6,0.6,0.6],...
                        'MarkerSize',5)

                    levelTrack = 0:size(plevel,1)-1;
                end
                plot(obj.sample(pList),...
                    zeros(length(pList),1),...
                    'o',...
                    'MarkerEdgeColor',[1,0,0],...
                    'MarkerFaceColor',[1,0,0],...
                    'MarkerSize',5)

                % boundries of sample markers
                plot(obj.sample(pdiff{1,1}(:,1)),...
                    (size(plevel,1)-1)*...
                    ones(size(pdiff{1,1}(:,1),1),1),...
                    'o',...
                    'MarkerEdgeColor',[0,0,0],...
                    'MarkerFaceColor',[0.6,0.6,0.6],...
                    'MarkerSize',8)
                % new partion markers
                for k = 2:size(pdiff,1)
                    plot(obj.sample(pdiff{k,1}(:,1)),...
                        (size(plevel,1)-k)*...
                        ones(size(pdiff{k,1}(:,1),1),1),...
                        'o',...
                        'MarkerEdgeColor',[0,0,0],...
                        'MarkerFaceColor',[0,0,0],...
                        'MarkerSize',8)
                end

                str = cell(1,size(levelTrack,2));
                for ii = 1:length(levelTrack)
                    str{ii} = sprintf('%1.0f',levelTrack(end+1-ii));
                end

                yticks(levelTrack)
                yticklabels(str)
                xlabel('x Range','Interpreter','latex')
                ylabel('Tree Level','Interpreter','latex')
                xlim([x_min x_max])

                bp = gca;
                if obj.save_figs
%                     saveas(bp, fullfile(fig_dir, [fig_name, '.png']))
                    saveas(bp, fullfile(fig_dir, [fig_name, '.fig']))
                    exportgraphics(bp,fullfile(fig_dir, [fig_name, '.png']),'Resolution',res)
                end
            end


        end

        function [brL,brR,partition] = gss(obj, sample)

            Ns = length(sample);
            bLeft = 1 + obj.binMin;
            bRight = Ns - obj.binMin;
            goldenR = (1+sqrt(5))/2;
            partition = ceil(Ns/2);

            while bRight-bLeft > 2

                % define partitions ---------------------------------------
                bLeft1 = partition;
                bRight1 = Ns - partition;
                rRight1 = obj.get_ratio(obj, sample(partition:end));
                rLeft1 = obj.get_ratio(obj, sample(1:partition));
                dxbrC = abs(obj.br_product(obj, bLeft1, rLeft1, Ns)-obj.br_product(obj, bRight1, rRight1, Ns));
                % --------------------
                leftPar = partition - 1;
                bLeft2 = leftPar;
                bRight2 = Ns - leftPar ;
                rRight2 = obj.get_ratio(obj, sample(leftPar:end));
                rLeft2 = obj.get_ratio(obj, sample(1:leftPar));
                dxbrL = abs(obj.br_product(obj, bLeft2,rLeft2,Ns)-obj.br_product(obj, bRight2,rRight2,Ns));
                % --------------------
                rightPar = partition + 1;
                bLeft3 = rightPar;
                bRight3 = Ns - rightPar;
                rRight3 = obj.get_ratio(obj, sample(rightPar:end));
                rLeft3 = obj.get_ratio(obj, sample(1:rightPar));
                dxbrR = abs(obj.br_product(obj, bLeft3,rLeft3,Ns)-obj.br_product(obj, bRight3,rRight3,Ns));
                % --------------------
                dxCR = dxbrC - dxbrR;
                dxLC = dxbrL - dxbrC;
                dxLR = dxbrL - dxbrR;
                
                % dxbrC = dxbrL or dxbrR = dxbrC or dxbrC is smallest
                if dxLC == 0 || dxCR == 0 || (dxLC > 0 && dxCR < 0) || bRight-bLeft <= 2
                    brL = obj.br_product(obj, bLeft1,rLeft1,Ns);
                    brR = obj.br_product(obj, bRight1, rRight1,Ns);
                    partition = bLeft1;
                    return
                end

                % parition update conditions ------------------------------
                 % dxbrR is smallest
                if dxCR >= 0 && dxLR > 0
                    bLeft = leftPar;
                    % Shrink --->
                    partition = round((bLeft+bRight*goldenR)/(1+goldenR));
                % dxbrL is smallest
                elseif dxLC <= 0 && dxLR < 0
                    bRight = rightPar;
                    % <--- Shrink
                    partition = round((bLeft*goldenR+bRight)/(1+goldenR));
                end
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

                % define partitions ---------------------------------------
                bLeft1 = partition;
                bRight1 = Ns - partition;
                rRight1 = obj.get_ratio(obj, sample(partition:end));
                rLeft1 = obj.get_ratio(obj, sample(1:partition));
                % --------------------
                leftPar = partition - 1;
                bLeft2 = leftPar;
                bRight2 = Ns - leftPar ;
                rRight2 = obj.get_ratio(obj, sample(leftPar:end));
                rLeft2 = obj.get_ratio(obj, sample(1:leftPar));
                % --------------------
                rightPar = partition + 1;
                bLeft3 = rightPar;
                bRight3 = Ns - rightPar;
                rRight3 = obj.get_ratio(obj, sample(rightPar:end));
                rLeft3 = obj.get_ratio(obj, sample(1:rightPar));
                % --------------------

                % MINIMIZE ABS()
                dxbrC = abs(obj.br_product(obj, bLeft1, rLeft1, Ns)-obj.br_product(obj, bRight1, rRight1, Ns));
                dxbrL = abs(obj.br_product(obj, bLeft2,rLeft2,Ns)-obj.br_product(obj, bRight2,rRight2,Ns));
                dxbrR = abs(obj.br_product(obj, bLeft3,rLeft3,Ns)-obj.br_product(obj, bRight3,rRight3,Ns));

                % MINIMIZE MAX()
%                 dxbrC = max(obj.br_product(obj, bLeft1, rLeft1, Ns), obj.br_product(obj, bRight1, rRight1, Ns));
%                 dxbrL = max(obj.br_product(obj, bLeft2,rLeft2,Ns), obj.br_product(obj, bRight2,rRight2,Ns));
%                 dxbrR = max(obj.br_product(obj, bLeft3,rLeft3,Ns),obj.br_product(obj, bRight3,rRight3,Ns));

                % MINIMIZE XI^2
%                 dxbrC = (obj.br_product(obj, bLeft1, rLeft1, Ns))^2+(obj.br_product(obj, bRight1, rRight1, Ns))^2;
%                 dxbrL = (obj.br_product(obj, bLeft2,rLeft2,Ns))^2+(obj.br_product(obj, bRight2,rRight2,Ns))^2;
%                 dxbrR = (obj.br_product(obj, bLeft3,rLeft3,Ns))^2+(obj.br_product(obj, bRight3,rRight3,Ns))^2;

                dxCR = dxbrC - dxbrR;
                dxLC = dxbrL - dxbrC;
                dxLR = dxbrL - dxbrR;

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
            % change one binary 1 value to 0 to create mask for difference array
            % indexes where mask is true i.e. 1
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

