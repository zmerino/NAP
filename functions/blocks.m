classdef blocks < NSE % inherit NSE properties i.e. p vector and max block size
    properties (Constant)
        % plot anything
        plt_figs = false;
        % plot xi values per parition or number of branches per level
        plt_xi = false;
        plt_tree = true;
    end
    properties
        sample;
        Ns;
        window;
        binMin;
        binNs;
    end
    methods(Static)
        function [j,nBlocks,kBlockLower,kBlockUpper,kList,T,BRlevel,BR0] = ...
                bin_width_size(obj)

            obj.Ns = length(obj.sample);

            if( obj.Ns < obj.binNs )
                nBlocks = 1;
            else
                % determine number of blocks and their sizes
                [kList,T,BRlevel,BR0] = obj.r_tree(obj);
                % remove tiny left and right
                if( kList(1) < obj.binNs )
                    kList = kList(2:end);
                end
                if( obj.Ns - kList(end) + 1 < obj.binNs )
                    kList = kList(1:end-1);
                end
                j = length(kList);
                nBlocks = 2*j + 1;
                kBlockLower = zeros(1,nBlocks);
                kBlockUpper = zeros(1,nBlocks);
            end
        end

        function [pList,T,BRlevel,BR0] = r_tree(obj)
            % parameters vector--------------------------------------------
            % p = [b, r, n, n-for-T, T-scale, window, maximum levels];

            % create T threshold-------------------------------------------
            %             xi = ( B^p(1)*R^p(2) ) / Ns^p(3);
            T = obj.p(5)*obj.Ns^obj.p(4);
            obj.window = ceil(obj.p(6)*obj.Ns^obj.p(7));
            maxLevel = obj.p(8);

            % minimum blocksize--------------------------------------------
            obj.binMin = ceil(2*obj.window);
            if obj.binMin > obj.Ns
                error('sample size too small for window size')
            end

            % track number of branches per level
            nbranch = 1;
            % set end points of sample length as partition left (pL) and
            % partion right (pR)
            pL = 1;
            pR = obj.Ns;
            % initialize vector to track all created partitions
            pList = [pL pR];
            % track every attempted partition for all levels
            plevel = {{[1;obj.Ns]}};
            % initialize array to track newly created partitions
            % for plotting purposes
            pdiff = {[1;obj.Ns]};
            % calcualte inital BR of intire sample
            B0 = obj.Ns;
            R0 = obj.get_ratio(obj, obj.sample);
            BR0 = obj.br_product(obj, B0, R0, obj.Ns);

            % clear array to hold all BR values per block per level
            BRlevel = {BR0};

            % beggin level loop
            %             if BR0 < T
            if BR0 > T
                for jj = 1:maxLevel
                    % vector to hold all attempted partitions per level
                    plevHold = [];
                    % vector to hold all BR values per level
                    BRHold = [];
                    % beggin branch loop
                    for b = 1:nbranch
                        % define block size (B)
                        B = pList(b+1) - pList(b);
                        if obj.binMin >= B
                            continue
                        end
                        % update left
                        boundryL = obj.binMin + 1;
                        boundryR = B - obj.binMin;
                        % block to small for for minimization
                        % given window size
                        if boundryR - boundryL < 3
                            break;
                        end

                        % golden ration bifraction minimization
                        [brL,brR,partition] = obj.min_br_gold(obj, obj.sample(pList(b):pList(b+1)));

                        % plot br and r for various paritions of whole
                        % sample
                        if jj == 1 && obj.plt_figs && obj.plt_xi
                            r_left = [];
                            r_right = [];
                            b_left = [];
                            b_right = [];
                            br_left = [];
                            br_right = [];
                            x_partition_left = [];
                            for n = 2*obj.window:Ns-2*obj.window
                                B0 = n;
                                b_left = [b_left, B0];
                                R0 = obj.get_ratio(sample(1:B0),obj.window);
                                r_left = [r_left, R0];
                                BR0 = obj.br_product(obj, B0, R0, obj.Ns);
                                br_left = [br_left, BR0];
                                x_partition_left = [x_partition_left, n];
                            end
                            x_partition_right = [];
                            for n = 2*obj.window:Ns-2*obj.window
                                B0 = n;
                                b_right = [b_right, B0];
                                R0 = obj.get_ratio(sample(end-B0:end),obj.window);
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

                            figure('Name', 'BR verse partition')
                            hold on;
                            plot(x_partition_left, br_left, '--b')
                            plot(x_partition_right, br_right, '--r')
                            plot(x_partition_right, abs(br_right-flip(br_left)), '-k')
                            xlabel('Partition', Interpreter='latex')
                            ylabel('$| \delta \xi |$', Interpreter='latex')
                        end

                        % find minimum BR for two newly created blocks
                        BR = min(brL,brR);
                        % update block boundaries of sample with correctly
                        % placed partition
                        newPar = pList(b) + partition;
                        % hold all BRs per level for later evaluation
                        BRHold = [BRHold,BR];
                        % un-balanced tree-----------------------------
                        %                         if BR <= T
                        if BR >= T
                            plevHold = [plevHold, newPar];
                        end
                    end
                    pList = sort([pList plevHold]);
                    % update nbranch
                    nbranch = length(pList) - 1;
                    % assign partition list to array for plotting
                    plevel{jj+1,1} = {pList'};
                    % exit for special cases where B < binMin
                    % or bRight - bLeft < 3
                    if isempty(BRHold)
                        break;
                    end
                    % assign BR per level to array for plotting
                    BRlevel{jj+1,1} = BRHold;
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
                pList = pList';

                %SPLITTING ROUTINE FOR LARGE SAMPLES ======================
                % vector to to add new partitons too
                LargNcheck = pList;
                % vector to hold updated partition list
                holder = pList;
                % while loop flag
                runSplit = true;

                while runSplit
                    % triggers exit flag for while loop
                    splitCount = 0;
                    % loop over modified partition list (holder)
                    for k = 1:length(holder)-1
                        % calcualte difference
                        diff = holder(k+1)-holder(k);
                        % add partion between elements when diff > max_bs
                        if holder(k+1)-holder(k) > obj.max_bs
                            split = floor((holder(k+1)-holder(k))/2);
                            % update new partiton list
                            LargNcheck = [LargNcheck,LargNcheck(k)+ split];
                            % update counter: number of found splits
                            splitCount = splitCount + 1;
                        end
                    end
                    LargNcheck = sort(LargNcheck);
                    holder = LargNcheck;
                    % if no splits exit routine
                    if splitCount == 0
                        runSplit = false;
                    end
                end
                pList = LargNcheck;
            end

            % FIGURES -----------------------------------------------------
            if obj.plt_figs && obj.plt_tree
                figure('Name','br values per level 1')
                hold on
                plot(0:size(plevel,1)-1,log(T*ones(size(plevel,1))), '-r');
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
                ylabel('ln(BR)')
                xlabel('Tree Level')
                legend('Threshold')

                %                 figure('Name','br values per level')
                %                 hold on
                %                 plot(0:size(plevel,1)-1,T*ones(size(plevel,1)), '-r');
                %                 for k = 2:size(BRlevel,1)
                %                     plot((k-2)*ones(size(BRlevel{k,1}(1,:),1),1),...
                %                         BRlevel{k,1}(1,:),...
                %                         'o',...
                %                         'MarkerEdgeColor',[0,0,0],...
                %                         'MarkerFaceColor',[0,0,0],...
                %                         'MarkerSize',4,'DisplayName','none')
                %                     levelTrack = 1:size(plevel,1);
                %                 end
                %
                %                 str = cell(1,size(levelTrack,2));
                %                 for ii = 1:length(levelTrack)-1
                %                     str{ii} = sprintf('%1.0f',levelTrack(ii));
                %                 end
                %                 xticks(levelTrack)
                %                 xticklabels(str)
                %                 ylabel('BR')
                %                 xlabel('Tree Level')
                %                 legend('Threshold')

                figure('Name','tree branching')
                subplot(2,1,1)
                histogram(obj.sample)
                ylabel('Number of Data Points','Interpreter','latex')

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
                % final partition markers
                % plot(obj.sample(plevel{end,1}{1,1}(:,1)),...
                %     zeros(size(plevel{end,1}{1,1}(:,1),1),1),...
                %     'o',...
                %     'MarkerEdgeColor',[1,0,0],...
                %     'MarkerFaceColor',[1,0,0],...
                %     'MarkerSize',5)
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
            end


        end

        function [brL,brR,partition] = min_br_gold(obj, sample)

            % function definition
            %             br = @br_product;
            %             sample = sort(sample);
            Ns = length(sample);
            bLeft = 1 + obj.binMin;
            bRight = Ns - obj.binMin;
            goldenR = (1+sqrt(5))/2;
            dxbr = [];
            loopCount = 0;
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

                % parition update conditions ------------------------------
                % dxbrR is smallest
                if dxCR >= 0 && dxLR > 0
                    bLeft = leftPar;
                    if bRight-bLeft < 2
                        dxbr = horzcat(dxbr,[dxbrC;bLeft1;rLeft1;bRight1;rRight1]);
                    else
                        % Shrink --->
                        dxbr = horzcat(dxbr,[dxbrR;bLeft3;rLeft3;bRight3;rRight3]);
                        partition = round((bLeft+bRight*goldenR)/(1+goldenR));
                    end
                end
                % dxbrL is smallest
                if dxLC <= 0 && dxLR < 0
                    bRight = rightPar;
                    if bRight-bLeft < 2
                        dxbr = horzcat(dxbr,[dxbrC;bLeft1;rLeft1;bRight1;rRight1]);
                    else
                        % <--- Shrink
                        dxbr = horzcat(dxbr,[dxbrL;bLeft2;rLeft2;bRight2;rRight2]);
                        partition = round((bLeft*goldenR+bRight)/(1+goldenR));
                    end
                end
                % dxbrC is smallest
                if dxLC > 0 && dxCR < 0
                    dxbr = horzcat(dxbr,[dxbrC;bLeft1;rLeft1;bRight1;rRight1]);
                    break
                end
                % dxbrC = dxbrL = dxbrR = 0
                if dxLC == 0 || dxCR == 0 || dxLR == 0
                    dxbr = horzcat(dxbr,[dxbrC;bLeft1;rLeft1;bRight1;rRight1]);
                    break
                end
                loopCount = loopCount + 1;
            end

            dxbr = dxbr';
            brL = obj.br_product(obj, dxbr(end,2),dxbr(end,3),Ns);
            brR = obj.br_product(obj, dxbr(end,4),dxbr(end,5),Ns);
            partition = dxbr(end,2);
        end

        function r = get_ratio(obj, sample)
            n = length(sample);
            dx = zeros(1,n-1);
            dx(1:n-1) = sample(2:n) - sample(1:n-1);
            % REQUIRED to sort the magnitude of differences
            dx = sort(dx');
            dxMin = mean(dx(1:obj.window));
            dxMax = mean(dx(end-obj.window:end));
            % optional ratio condition
            % -------------------------------------------------------------
            % r -> 0, dxMin << dxMax and r -> 1, dxMin ~ dxMax
            % -------------------------------------------------------------
            %             r = dxMin/dxMax;
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

