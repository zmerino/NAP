function [pdfPoints, pdfEst] = PDFe(r)
    nVariables = size(r, 2);
    nSamples = size(r, 1);
    if nVariables == 1
        parm.debug = 0;
        [~, pdfPoints, pdfEst] = EstimatePDF(r, parm);
        pdfPoints = pdfPoints';
        return
    end
    maxGrid = 100;
    minGrid = 3;
    minSamples = 100;

    nGrids = ceil((nSamples / minSamples) ^ (1 / (nVariables - 1)));
    if nGrids > maxGrid
        nGrids = maxGrid;
    elseif nGrids < minGrid
        nGrids = minGrid;
    end
  
    [jp, x] = EstimatePDFmv(r, nSamples, nVariables, nGrids);
    
    if nVariables == 2
        pdfEst = reshape(jp, nGrids, nGrids);
    elseif nVariables == 3
        pdfEst = reshape(jp, nGrids, nGrids, nGrids);
    elseif nVariables == 4
        pdfEst = reshape(jp, nGrids, nGrids, nGrids, nGrids);
    elseif nVariables == 5
        pdfEst = reshape(jp, nGrids, nGrids, nGrids, nGrids, nGrids);
    end
    pdfPoints = transpose(reshape(x, nGrids, nVariables));
    pdfEst = permute(pdfEst, flip(1:nVariables));    

end