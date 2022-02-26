function [dxbr1,brL,brR,partition] = minimizeBRgold(sample,window,p,binMin)

% function definition
br = @brProduct;

% sample size
sample = sort(sample);
Ns = length(sample);

% disp('Ns')
% disp(Ns)


bLeft = 1 + binMin;
bRight = Ns - binMin;

goldenR = (1+sqrt(5))/2;
jiggle = 1;
dxbr = [];
loopCount = 0;

tic
partition = ceil(Ns/2);

while bRight-bLeft > 2
   

    bLeft1 = partition ;
    bRight1 = Ns - partition ;
    
    rRight1 = getRation(sample(partition:end),window);
    rLeft1 = getRation(sample(1:partition),window);
    
    dxbrC = abs(br(bLeft1,rLeft1,p,Ns)-br(bRight1,rRight1,p,Ns));
    
    % ----------------------------------------------------------
    leftPar = partition - jiggle;
    
    bLeft2 = leftPar;
    bRight2 = Ns - leftPar ;
    
    rRight2 = getRation(sample(leftPar:end),window);
    rLeft2 = getRation(sample(1:leftPar),window);
    
    dxbrL = abs(br(bLeft2,rLeft2,p,Ns)-br(bRight2,rRight2,p,Ns));
    
    % ----------------------------------------------------------
    rightPar = partition + jiggle;
    
    bLeft3 = rightPar;
    bRight3 = Ns - rightPar ;
    rRight3 = getRation(sample(rightPar:end),window);
    rLeft3 = getRation(sample(1:rightPar),window);
    
    dxbrR = abs(br(bLeft3,rLeft3,p,Ns)-br(bRight3,rRight3,p,Ns));
    
    % ----------------------------------------------------------
    dxCR = dxbrC - dxbrR;
    dxCL = dxbrC - dxbrL;
    dxLR = dxbrL - dxbrR;
    
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
    if dxCL >= 0 && dxLR < 0
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
    if dxCL < 0 && dxCR < 0
        dxbrFinal = dxbrC;
        dxbr = horzcat(dxbr,[dxbrC;bLeft1;rLeft1;bRight1;rRight1]);
        break
    end
    
    loopCount = loopCount + 1;
end
timeRandom = toc;

% disp(' ')
% disp('Program Results =======================================')
dxbr = dxbr';
% dxbrDisplay1 = [dxbr(:,1),dxbr(:,2)];
dxbr1 = dxbr(end,1);
bL1 = partition;


% dxbr(end,2),dxbr(end,3)
% dxbr(end,4),dxbr(end,5)

% brL = br(bLeft1,rLeft1,p,Ns);
% brR = br(bRight1,rRight1,p,Ns);
brL = br(dxbr(end,2),dxbr(end,3),p,Ns);
brR = br(dxbr(end,4),dxbr(end,5),p,Ns);

partition = dxbr(end,2);
bR1 = Ns - bL1;

end







