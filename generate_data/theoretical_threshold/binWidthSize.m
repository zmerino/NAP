function [T0,j,nBlocks,kBlockLower,kBlockUpper,kList] = binWidthSize(Ns,binNs,sample,filename,savePNG,lowLim,upLim,p)
if( Ns < binNs )                                        % => use one block
    nBlocks = 1;
else
    % determine number of blocks and their sizes
    [kList,T0] = Rtree(sample,lowLim,upLim,savePNG,filename,p);
    %kList = sort(kList);
    % ------------------------------------- remove tiny left and right
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