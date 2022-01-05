classdef block_definition
    methods(Static)
        function [j,nBlocks,kBlockLower,kBlockUpper,kList,T,BRlevel,BR0] = ...
                bin_width_size(Ns,binNs,sample,filename,lowLim,upLim,p)
            if( Ns < binNs )
                nBlocks = 1;
            else
                % determine number of blocks and their sizes
                [kList,T,BRlevel,BR0] = block_definition.r_tree(sample,lowLim,upLim,filename,p);
                % remove tiny left and right
                if( kList(1) < binNs )
                    kList = kList(2:end);
                end
                if( Ns - kList(end) + 1 < binNs )
                    kList = kList(1:end-1);
                end
                j = length(kList);
                nBlocks = 2*j + 1;
                kBlockLower = zeros(1,nBlocks);
                kBlockUpper = zeros(1,nBlocks);
            end
        end
        
        function [pList,T,BRlevel,BR0] = r_tree(sample,lowLim,upLim,filename,p)
            debug = false;  %<- true/false display debug yes/no
            % Initialize variables
            % create functions
            %             br = @br_product;
            %             r = @get_ratio;
            % length of sample
            Ns = length(sample);
            % parameters vector--------------------------------------------
            % p = [b, r, n, n-for-T, T-scale, window, maximum levels];
            % BR product function parameters-------------------------------
            % br = @(b,r,p,n) (b^p(1)*r^p(2))/n^p(3);
            %--------------------------------------------------------------
            % % p = [1,0.5,1,0.25,8,ceil(0.125*Ns^0.5),20];
            % % % p original
            % % p = [1,0.5,1,0.25,8,ceil(0.0625*Ns^0.5),40];
            % % % p realy low T; controled by max level
            % % p = [1,0.5,1,0.33,2,ceil(0.0625*Ns^0.5),40];
            % % % % p = [1,0.5,1,0.25,8,ceil(0.0125*Ns),7];
            % p = [1,1,1,0.25,8,ceil(0.125*Ns^0.5),20];
            p = [1,0.55,1,0.33,2,ceil(0.0625*Ns^0.5),40];
            p = [1,0.55,1,0.33,0.0002,6*ceil(Ns^0.5),40];
            p = [1,0.55,1,0.33,0.5,6*ceil(Ns^0.5),40];
            p = [0,1,0,0.33,0.5,6*ceil(Ns^0.5),40];
            % create T threshold-------------------------------------------
            %T = 0.07%p(5)*Ns^p(4)
% %             a = 1;
% %             b = 0.000009999999999999999;
% %             c = 27300;
% %             d = -.3;%0.173;
            a = 1;
            b = 0.000012;
            c = -78000;
            d = 0;
            T = (1/(a+exp(b*(Ns + c)))) + d;
            T = 2.6741*Ns^(-0.626) + 0.003;
            T = 2.6741*Ns^(-0.626) + 0.005;
            % window parameter for top/bottom points to average------------
            window = p(6);
            
            %disp(['Window percentage: ',num2str((window/Ns)*100)])
            
%             pause
            
            % maximum number of potential levels/splits--------------------
            maxLevel = p(7);
            % minimum blocksize--------------------------------------------
            binMin = ceil(2*window);
            if binMin > Ns
                error('sample size too small for window size')
            end
            
            % track number of branches per level
            nbranch = 1;
            % set end points of sample length as partition left (pL) and
            % partion right (pR)
            pL = 1;
            pR = Ns;
            % initialize vector to track all created partitions
            pList = [pL pR];
            % track every attempted partition for all levels
            plevel = {{[1;Ns]}};
            % initialize array to track newly created partitions
            % for plotting purposes
            pdiff = {[1;Ns]};
            % calcualte inital BR of intire sample
            B0 = Ns;
            
            R0 = block_definition.get_ratio(sample,window);
            BR0 = block_definition.br_product(B0,R0,p,Ns);
            
%             pause(1)

            %{
            %             T = p(5)*BR0*Ns^p(4);
            %             T = p(5)*(BR0*Ns)^p(4);
            %             T = 2*(BR0/Ns)^(0.5); % <- way too small
            %             T = 40*BR0*Ns^(-1); % <- works descently: T a little to large for small n and grows too slowly
            %             T = 600*Ns^(1) + 2^(18); % <- very good
            %             T = 4000*Ns^(1) + 2^(17);
            %             T = 8000*Ns^(1) + 2^(17);
            %             T = 2^(-8)*(22700*Ns^(1.25)-40000); %<- good but issues with failures for stable/pareto
            %
            % %             a = 0.44;
            % %             b = 1.35;
            % %             c = 16;
            % %             d = 6.8;
            % good but causes issue with block overlap for pareto and fails
            % for mixed-unirom
            %             a = 0.44;
            %             b = 1.35;
            %             c = 16;
            %             d = 9;
            
            %             a = 0.44;
            %             b = 1.5;
            %             c = 100;
            %             d = 9;
            %
            %             a = 0.104;
            %             b = 1.9;
            %             c = 100;
            %             d = 9;
            % %
            %             a = 0.44;
            %             b = 1.4;
            %             c = 100;
            %             d = 6.3;
            
            %             a = 0.44;
            %             b = 1.3;
            %             c = 100;
            %             d = 5;
            
            %             T = 2^(d)*(a*Ns^(b)+c);
            
            %-----------
            %             F = 1;
            %             G = 1;
            %             H = 0.7;
            %             I = -12.3;
            %
            %             e = F/(G+exp(-H*(Ns+I)));
            %             ep = 1 - e;
            %
            %             a_large = 0.351;
            %             b_large = 1.8;
            %             c_large = 100;
            %             d_large = 9;
            %
            %             a_small = 0.44;
            %             b_small = 1.35;
            %             c_small = 100;
            %             d_small = 9;
            %
            %             T_small = 2^(d_small)*(a_small*Ns^(b_small)+c_small);
            %             T_large = 2^(d_large)*(a_large*Ns^(b_large)+c_large);
            %
            %             T = T_small*e + T_large*ep;
            %---------------
            %             T = 16000*Ns^(1) + 2^(13);
            
            %             T = 1000*BR0*Ns^(-1.5); % <- went wrong direction from ^^
            %             T = 1*BR0*Ns^(-0.5);
            %             T = 20*BR0*Ns^(-0.75);
            %             T = 4*Ns^(1) + 100000;
            %             T = 0.01*BR0;
            
            %              disp(['Ns: ', num2str(Ns)])
            %              disp(['R0: ', num2str(R0)])
            %              disp(['B0: ', num2str(B0)])
            %              disp(['BR0: ', num2str(BR0)])
            %              disp(['BR0Ns: ', num2str(BR0*Ns^(-1))])
            %              disp(['window: ', num2str(window)])
            %              disp(['binMin: ', num2str(binMin)])
            %              disp(['T: ', num2str(T)])
            %              disp(' ')
            %}
            
            % display useful variable values
            if debug
                disp(['Ns: ', num2str(Ns)])
                disp(['R0: ', num2str(R0)])
                disp(['B0: ', num2str(B0)])
                disp(['BR0: ', num2str(BR0)])
                disp(['window: ', num2str(window)])
                disp(['binMin: ', num2str(binMin)])
                disp(['T: ', num2str(T)])
                disp(' ')
            end
            
            % clear array to hold all BR values per block per level
            BRlevel = {BR0};
            
            % beggin level loop
            if BR0 < T
                for jj = 1:maxLevel
                    if debug
                        disp(' ')
                        disp(['\\\\\\\\\\\\\\\ START LEVEL: ', num2str(jj)])
                    end
                    % vector to hold all attempted partitions per level
                    plevHold = [];
                    % vector to hold all BR values per level
                    BRHold = [];
                    % beggin branch loop
                    for b = 1:nbranch
                        if debug
                            disp(['-------------- START branch: ', num2str(b)])
                        end
                        % define block size (B)
                        B = pList(b+1) - pList(b);
                        if binMin >= B
                            if debug
                                disp(['block size (', num2str(B),...
                                    ') smaller than: ', num2str(binMin)])
                                disp(['End this branch'])
                            end
                            continue
                        end
                        % update left
                        boundryL = binMin + 1;
                        boundryR = B - binMin;
                        % block to small for for minimization
                        % given window size
                        if boundryR - boundryL < 3
                            break;
                        end
                        
                        % plot r,b,br per partition
%                          block_definition.minimizeBRdiff(sample(pList(b):pList(b+1)),window,p,binMin,filename,jj);
                        
                        % golden ration bifraction minimization
                        [dxBR,brL,brR,partition] = block_definition.min_br_gold(sample(pList(b):pList(b+1)),window,p,binMin,filename);
                        
                       
                        % find minimum BR for two newly created blocks
                        BR = min(brL,brR);
                        % update block boundaries of sample with correctly
                        % placed partition
                        newPar = pList(b) + partition;
                        % hold all BRs per level for later evaluation
                        BRHold = [BRHold,BR];
                        if debug
                            disp(['BR: ',num2str(BR)])
                            disp(['T:  ', num2str(T)])
                        end
                        % un-balanced tree-----------------------------
                        if BR <= T
                            plevHold = [plevHold, newPar];
                            if debug
                                disp(['partition1: ',num2str(newPar)])
                            end
                        end
                    end
                    if debug
                        disp(['BR: ',num2str(BR)])
                        disp(['T:  ', num2str(T)])
                    end
                    pList = [pList plevHold];
                    pList = sort(pList);
                    % update nbranch
                    nbranch = length(pList) - 1;
                    % assign partition list to array for plotting
                    plevel{jj+1,1} = {pList'};
                    % exit for special cases where B < binmin
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
                if debug
                    disp(' ')
                    disp('********* pList FINAL ANSWER')
                    disp(pList)
                    %                     disp(['Elapse time: ',num2str(endTime),'s'])
                    disp('******************************')
                end
                pList = pList';
                sample = sort(sample);
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
                        % add partion between elements when diff>20000
                        if holder(k+1)-holder(k) > 100000
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
                endTime = toc;
                pList = LargNcheck;
                if debug
                    disp('****************** pList SPLIT')
                    disp(pList')
                    disp(['Elapse time: ',num2str(endTime),'s'])
                    disp('******************************')
                end
                %----------------------------------------------------------
            end
           
        end
        
        function [dxbr1,brL,brR,partition] = min_br_gold(sample,window,p,binMin,filename)
            
            % function definition
            %             br = @br_product;
            sample = sort(sample);
            Ns = length(sample);
            bLeft = 1 + binMin;
            bRight = Ns - binMin;
            goldenR = (1+sqrt(5))/2;
            jiggle = 1;
            dxbr = [];
            loopCount = 0;
            partition = ceil(Ns/2);
            
%             disp('enter while loop +++++++++++++++++++++++++++++++++++++++++++++++++++++++')
            testTrip = 100;
            testx = zeros(1,testTrip);
            testy = zeros(1,testTrip);
            
            while bRight-bLeft > 2
                
%                 if  loopCount >= testTrip
%                     dxbr = horzcat(dxbr,[dxbrC;bLeft1;rLeft1;bRight1;rRight1]);
%                     break;
%                 end
                
                bLeft1 = partition;
                bRight1 = Ns - partition; 
                rRight1 = block_definition.get_ratio(sample(partition:end),window);
                rLeft1 = block_definition.get_ratio(sample(1:partition),window);
                dxbrC = abs(block_definition.br_product(bLeft1,rLeft1,p,Ns)-block_definition.br_product(bRight1,rRight1,p,Ns));
                % ----------------------------------------------------------
                leftPar = partition - jiggle;
                bLeft2 = leftPar;
                bRight2 = Ns - leftPar ;
                rRight2 = block_definition.get_ratio(sample(leftPar:end),window);
                rLeft2 = block_definition.get_ratio(sample(1:leftPar),window);
                dxbrL = abs(block_definition.br_product(bLeft2,rLeft2,p,Ns)-block_definition.br_product(bRight2,rRight2,p,Ns));
                % ----------------------------------------------------------
                rightPar = partition + jiggle;
                bLeft3 = rightPar;
                bRight3 = Ns - rightPar; 
                rRight3 = block_definition.get_ratio(sample(rightPar:end),window);
                rLeft3 = block_definition.get_ratio(sample(1:rightPar),window);
                dxbrR = abs(block_definition.br_product(bLeft3,rLeft3,p,Ns)-block_definition.br_product(bRight3,rRight3,p,Ns));
                % ----------------------------------------------------------
                dxCR = dxbrC - dxbrR;
                dxLC = dxbrL - dxbrC;
                dxLR = dxbrL - dxbrR;
%                 disp('--------------------------------------------------')
                
                % dxbrR is smallest
                if dxCR >= 0 && dxLR > 0
                    bLeft = leftPar;
                    %--------------------------------------------------
                    if bRight-bLeft < 2
                        dxbr = horzcat(dxbr,[dxbrC;bLeft1;rLeft1;bRight1;rRight1]);
                    else
                        dxbr = horzcat(dxbr,[dxbrR;bLeft3;rLeft3;bRight3;rRight3]);
                        partition = round((bLeft+bRight*goldenR)/(1+goldenR));
                        %             disp('Shrink --->')
                    end
                end
                
                % dxbrL is smallest
                if dxLC <= 0 && dxLR < 0
                    bRight = rightPar;
                    %--------------------------------------------------
                    if bRight-bLeft < 2
                        dxbr = horzcat(dxbr,[dxbrC;bLeft1;rLeft1;bRight1;rRight1]);
                    else
                        dxbr = horzcat(dxbr,[dxbrL;bLeft2;rLeft2;bRight2;rRight2]);
                        partition = round((bLeft*goldenR+bRight)/(1+goldenR));
                        %             disp('<--- Shrink')
                    end
                end
                
                % dxbrC is smallest
                if dxLC > 0 && dxCR < 0 
                    dxbrFinal = dxbrC;
                    dxbr = horzcat(dxbr,[dxbrC;bLeft1;rLeft1;bRight1;rRight1]);
                    break
                end
                
                % dxbrC = dxbrL = dxbrR = 0
                if dxLC == 0 || dxCR == 0 || dxCR == 0
                    dxbrFinal = dxbrC;
                    dxbr = horzcat(dxbr,[dxbrC;bLeft1;rLeft1;bRight1;rRight1]);
                    break
                end
                
                loopCount = loopCount + 1;
                
                testx(loopCount) = loopCount;
                testy(loopCount) = partition;
                
            end
            
            
%             disp(['exit while loop ',num2str(loopCount)])
            
%             figure('Name','test')
%             plot(testx,testy,'.r')
             
             
%             pause
            
            
            
            dxbr = dxbr';
            % dxbrDisplay1 = [dxbr(:,1),dxbr(:,2)];
            dxbr1 = dxbr(end,1);
            bL1 = partition;
            % dxbr(end,2),dxbr(end,3)
            % dxbr(end,4),dxbr(end,5)
            % brL = block_definition.br_product(bLeft1,rLeft1,p,Ns);
            % brR = block_definition.br_product(bRight1,rRight1,p,Ns);
            brL = block_definition.br_product(dxbr(end,2),dxbr(end,3),p,Ns);
            brR = block_definition.br_product(dxbr(end,4),dxbr(end,5),p,Ns);
            partition = dxbr(end,2);
            bR1 = Ns - bL1;
        end
        
        function r = get_ratio(sample,window)
            n = length(sample);
            dx = zeros(1,n-1);
            dx(1:n-1) = sample(2:n) - sample(1:n-1);
            dx = sort(dx');
            dxMin = mean(dx(1:window));
            dxMax = mean(dx(end-window:end));
%             r = dxMax/dxMin;
%             r = (dxMax/dxMin)*(dxMin+dxMax)
            r = dxMin/dxMax;
            
        end
        
        function [dxMin,dxMax] = get_ratioTest(sample,window)
            n = length(sample);
            dx = zeros(1,n-1);
            dx(1:n-1) = sample(2:n) - sample(1:n-1);
            dx = sort(dx');
            dxMin = mean(dx(1:window));
            dxMax = mean(dx(end-window:end));
%             %-- max
            %- min
%             dx(1:window)
%             pause
%             dx(end-window:end)
%             pause
            %             r = dxMax/dxMin;
            %             r = (dxMax/dxMin)*(dxMin+dxMax)
%             r = dxMin/dxMax;
            
        end
        
        function BR = br_product(b,r,p,n)
            BR = (b^p(1)*r^p(2))/n^p(3);
        end
        
        function [dxbr,brL,brR,partition] = minimizeBRdiff(sample,window,p,binMin,filename,j)
         
            
            Ns = length(sample);
            sample = sort(sample);
         
            step = 1;
            % step = ceil(0.1*Ns);
            % bL = 0;
            dxbrTrack = [];
            dxLMinTrack = [];
            dxRMinTrack = [];
            dxLMaxTrack = [];
            dxRMaxTrack = [];
            rLTrack = [];
            rRTrack = [];
            bRTrack = [];
            bLTrack = [];
            brLTrack = [];
            brRTrack = [];
            BRleft = [];
            BRright = [];
            trig = -1;
            rtreeFlag = false;
            while trig < 0
                
                %     bL = bL + binMin + step;
                bL = binMin + step;
                bR = Ns - step - binMin;
                
                if bR < binMin + 2
                    break;
                    rtreeFlag = true;
                end
                
                rL = block_definition.get_ratio(sample(1:bL),window);
                rR = block_definition.get_ratio(sample(bL:end),window);
                
                [RdxMin,RdxMax] = block_definition.get_ratioTest(sample(1:bL),window);
                [LdxMin,LdxMax] = block_definition.get_ratioTest(sample(bL:end),window);
                
                dxLMinTrack = [dxLMinTrack,LdxMin];
                dxRMinTrack = [dxRMinTrack,RdxMin];
                dxLMaxTrack = [dxLMaxTrack,LdxMax];
                dxRMaxTrack = [dxRMaxTrack,RdxMax];
                
                rLTrack = [rLTrack, rL];
                rRTrack = [rRTrack,rR];
                                
                dxbr = abs(block_definition.br_product(bL,rL,p,Ns)-block_definition.br_product(bR,rR,p,Ns));
                
                BRleft = [BRleft,block_definition.br_product(bL,rL,p,Ns)];
                BRright = [BRright,block_definition.br_product(bR,rR,p,Ns)];
                
                
                brLTrack = [brLTrack, block_definition.br_product(bL,rL,p,Ns)];
                brRTrack = [brRTrack, block_definition.br_product(bR,rR,p,Ns)];
                
                
                bLTrack = [bLTrack, bL];
                bRTrack = [bRTrack,bR];
                
                
                dxbrTrack = [dxbrTrack,dxbr];

                step = step + 1;
            end
            
            if j == 1
                
                
            end
            
            [dxbr,bL]= min(dxbrTrack);
            
            brL = brLTrack(bL);
            brR = brRTrack(bL);
            
            partition = bLTrack(bL);
            
            dxbrTrack = dxbrTrack';
            
        end
        
    end
end

