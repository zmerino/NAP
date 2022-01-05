
function [pList,T] = Rtree(sample,lowLim,upLim,saveImage,filename,p)

treeType =              false; %<- true/false ballanced/unbalanced
debug =                 false;  %<- true/false display debug yes/no
%% Initialize variables
% create functions
br = @brProduct;
r = @getRation;
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
% create T threshold-------------------------------------------
T = p(5)*Ns^p(4);
% T = p(5)*Log(Ns+1)^p(4);
% T = p(5)*log(Ns+1);
% window parameter for top/bottom points to average------------
window = p(6);
% maximum number of potential levels/splits--------------------
maxLevel = p(7);
% minimum blocksize--------------------------------------------
binMin = ceil(2*window);

if binMin > Ns
    error('sample size too small for window size')
end

% display useful variable values
if debug
    disp(['Ns: ', num2str(Ns)])
    disp(['T: ', num2str(T)])
    disp(['window: ', num2str(window)])
    disp(['binMin: ', num2str(binMin)])
    disp(' ')
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

% initial exit flag
exit = false;

% calcualte inital BR of intire sample
B0 = Ns;
R0 = r(sample,window);
BR0 = br(B0,R0,p,Ns);
% clear array to hold all BR values per block per level
BRlevel = {BR0};

tic
t = cputime;
% beggin level loop
if BR0 > T
    
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
            %                         given window size
            if boundryR - boundryL < 3
                break;
            end
            
            % golden ration bifraction minimization
            [dxBR,brL,brR,partition] = minimizeBRgold(sample(pList(b):pList(b+1)),window,p,binMin);
            
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
            if BR >= T
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
        disp(['Elapse time: ',num2str(endTime),'s'])
        disp('******************************')
    end
    pList = pList';
    sample = sort(sample);
    
    %SPLITTING ROUTINE FOR LARGE SAMPLES-----------------------
    
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
% FIGURES -----------------------------------------------------
figure('Name','br values per level')
hold on
plot(0:size(plevel,1),log(T*ones(size(plevel,1)+1)), '-r');
for k = 1:size(BRlevel,1)
    plot((k-1)*ones(size(BRlevel{k,1}(1,:),1),1),...
        log(BRlevel{k,1}(1,:)),...
        'o',...
        'MarkerEdgeColor',[0,0,0],...
        'MarkerFaceColor',[0,0,0],...
        'MarkerSize',4,'DisplayName','none')
    levelTrack = 1:size(plevel,1);
end

str = cell(1,size(levelTrack,2));
for ii = 1:length(levelTrack)
    str{ii} = sprintf('%1.0f',levelTrack(ii));
end
xticks(levelTrack)
xticklabels(str)
ylabel('ln(BR)')
xlabel('Tree Level')
legend('Threshold')
if saveImage
    binFileName = ['BR_',char(filename)];
    pngfile = strcat(binFileName,'.png');
    saveas(gcf,pngfile)
    figfile = strcat(binFileName,'.fig');
    saveas(gcf,figfile)
end
%--------------------------------------------------------------------------
figure('Name','tree branching')
subplot(2,1,1)
histogram(sample)
ylabel('Number of Data Points','Interpreter','latex')

subplot(2,1,2)
hold on
% branching level track markers
for k = 1:size(plevel,1)-1
    plot(sample(plevel{k,1}{1,1}(:,1)),...
        (size(plevel,1)-k)*...
        ones(size(plevel{k,1}{1,1}(:,1),1),1),...
        'o',...
        'MarkerEdgeColor',[0.6,0.6,0.6],...
        'MarkerFaceColor',[0.6,0.6,0.6],...
        'MarkerSize',5)
    
    levelTrack = 0:size(plevel,1)-1;
end
% final partition markers
%             plot(sample(plevel{end,1}{1,1}(:,1)),...
%                 zeros(size(plevel{end,1}{1,1}(:,1),1),1),...
%                 'o',...
%                 'MarkerEdgeColor',[1,0,0],...
%                 'MarkerFaceColor',[1,0,0],...
%                 'MarkerSize',5)

plot(sample(pList),...
    zeros(length(pList),1),...
    'o',...
    'MarkerEdgeColor',[1,0,0],...
    'MarkerFaceColor',[1,0,0],...
    'MarkerSize',5)

% boundries of sample markers
plot(sample(pdiff{1,1}(:,1)),...
    (size(plevel,1)-1)*...
    ones(size(pdiff{1,1}(:,1),1),1),...
    'o',...
    'MarkerEdgeColor',[0,0,0],...
    'MarkerFaceColor',[0.6,0.6,0.6],...
    'MarkerSize',8)
% new partion markers
for k = 2:size(pdiff,1)
    plot(sample(pdiff{k,1}(:,1)),...
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
if saveImage
    binFileName = ['Tree_',char(filename)];
    pngfile = strcat(binFileName,'.png');
    saveas(gcf,pngfile)
    figfile = strcat(binFileName,'.fig');
    saveas(gcf,figfile)
end
end

