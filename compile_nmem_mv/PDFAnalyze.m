function [py, y, cy, sqr] = PDFAnalyze(sample, varargin)
%PDFAnalyze computes a probability density estimate for a one-dimensional
%data sample and produces optional plots for analysis 
%
% [F, XI] = PDFAnalyze(X) Computes the density estimate of data in
% sample X.  F contains the density estimate at points XI.  The number of
% points and the relative spacing is determined automatically from the
% features of the data sample
%
% [F, XI, CDF, SQR] = PDFAnalyze(...) also returns the cumulative density
% and the scaled quantile residual for each sample data point.
%
% PDFAnalyze(...) with no output arguments produces a plot of the density estimate.
%
% [...] = PDFAnalyze(..., 'param1', 'val1, 'param2', 'val2', ...) specifies
% parameter name/value pairs to control the density estimation.  Valid
% parameters are as follows:
%
%   Parameter               Value
%   'PlotType'              Produces any combination of three plot types:
%                           'pdf'       probability density function for the
%                           'sqr'       scaled quantile residual
%                           'combined'  pdf and sqr plotted on one figure
%
%                           Multiple plot types occur with multiple name/value pairs specified
%
%                           The equation for Scaled Quantile Residual (SQR) is given by SQR = sqrt(N+2)*(u - uniform-u) where N is the number of data samples. 
%                           SQR plots are very useful as a diagnostic measure because they are sample size invariant and have universal characteristics 
%                           independent of the true PDF. The SQR plot type plots the SQR for each data sample by position, highlighting in red those that fall 
%                           outside of the expected 98% threshold.
%
%
%   'EstimationType'        The default estimation method is PDFEstimate.
%                           Additional KDE methods are available:
%                           'kde1'      built-in MATLAB function ksdensity
%                           'kde2'      Zdravko Botev (2020). Kernel Density Estimator 
%                                       (https://www.mathworks.com/matlabcentral/fileexchange/14034-kernel-density-estimator), 
%                                       MATLAB Central File Exchange. Retrieved March 17, 2020.
%
%   'distribution'          A two column matrix, [F, XI], representing a
%                           distribution to plot on the same figure as the
%                           estimate for use with 'pdf' plot type.  Useful
%                           for comparison to a known distribution.
%
%
%   Example: Plot the estimate of random sample for a Normal distribution
%   along with the true Normal distribution:
%       data = randn(1000, 1);
%       x = min(data):0.1:max(data);
%       f = normpdf(x);
%       d = [x(:), f(:)];
%       PDFAnalyze(data, 'distribution', d);
%
%   Example: Plot the scaled quantile residual (SQR) for an estimate of the
%   Normal distribution, showing confidence thresholds and uncertainty
%   estimates:
%       PDFAnalyze(randn(10000, 1), 'PlotType', 'sqr');
%
% See also ksdensity
%


    plotPDF = 0;
    plotExact = 0;
    plotSQR = 0;
    plotBoth = 0;
    
    estimatePDFE = 1;
    estimateKDE1 = 0;
    estimateKDE2 = 0;
    
    failed = false;
    greyColor = [0.65, 0.65, 0.65];
    redColor = [[0.65,0.078,0.18]];
    blueColor = [45/255, 66/255, 133/255];
    nargs = length(varargin);
    for j=1:2:nargs
        if strcmpi(varargin{j}, 'plotType')
            FigureSettings();
            parameters.adaptiveDx = false;              
            if (strcmpi(varargin{j+1}, 'pdf'))
                plotPDF = 1;
            elseif (strcmpi(varargin{j+1}, 'sqr'))
                plotSQR = 1;
            elseif (strcmpi(varargin{j+1}, 'combined'))
                plotBoth = 1;
            else
                disp(['Unknown value for parameter plotType: ', ': ', varargin{j+1}]);
            end        
        elseif strcmpi(varargin{j}, 'estimationType')
            if strcmpi(varargin{j+1}, 'kde1')
                estimateKDE1 = 1;
                estimatePDFE = 0;
            elseif strcmpi(varargin{j+1}, 'kde2')
                estimateKDE2 = 1;
                estimatePDFE = 0;  
            else                
                disp(['Unknown value for parameter estimationType: ', ': ', varargin{j+1}]);              
            end
        elseif strcmpi(varargin{j}, 'distribution')
            plotExact = 1;
            distribution = varargin{j+1};
        else            
            disp(['Unknown parameter: ', varargin{j}]);
        end
    end    
    
    if ~plotPDF && ~plotSQR && ~plotBoth
        if nargout == 0
            plotPDF = 1;
        end
    end
       
    if estimatePDFE
        parameters.SURDtarget = 40;                 % {number > 0 and <= 100}
        parameters.SURDmin = 5;                     % {number > 0}
        parameters.SURDmax = 100;                   % {number > SURDmin and <= 100}
        parameters.LagrangeMin = 1;                 % {number >= 1}
        parameters.LagrangeMax = 200;               % {number > LagrangeMin}
        parameters.partition = 0;                   % {any number, zero for no partitioning}
        parameters.debug = false;                   % {true, false}
        parameters.scoreType = 'QZ';                % {'QZ', 'LL'}
        parameters.outlierCutoff = 7;               % {number >= 0, zero to keep all outliers}
        try
            [failed, y, py, cy, sqr] = EstimatePDF(sample, parameters);
        catch ME
            failed = true;
            y = 0;
            py = 0;
            cy = 0;
            sqr = 0;
            disp(ME.message);
        end
    else
        n = length(sample);       
        if estimateKDE1
           [cy, ~] = ksdensity(sample, 'Function', 'cdf');
           [py, y] = ksdensity(sample);
           if isnan(y)
               failed = true;
           end    
        end
        if estimateKDE2
            [~, py, y, cy] = kde(sample, 2^12, [0 1]);
        end
        position = (1:n) ./ n; 
        [~, idxUnique, ~] = unique(cy);
        rCalc = interp1(y(idxUnique), cy(idxUnique), sort(sample));  
        residualCalc = rCalc - position;
        sqr = residualCalc * sqrt(length(rCalc) + 2);
        sqr = sqr';
    end
    
    if ~failed
        idxLeft  = find(sample > min(y));
        idxRight = find(sample < max(y));
        data = sort(sample(intersect(idxLeft, idxRight)));
        
        if plotPDF  
            hPDF = figure;
            fPDF = axes('parent', hPDF);
            if plotExact
                hold on;
                plot(distribution(:, 1), distribution(:, 2), 'color', greyColor, 'linewidth', 3.0);
                plot(fPDF, y, py, 'color', blueColor, 'linewidth', 2.0);
                legend('exact', 'estimate');
            else
                plot(fPDF, y, py, 'color', blueColor);
            end
            ylabel('PDF');                 
            xlabel('Sample Range');
        end
        
        if plotSQR  
            hSQR = figure;
            fSQR = axes('parent', hSQR);
            [topThreshold, bottomThreshold] = PlotBeta(data, false);
            n = length(sqr);
            dx = 1 / (n + 1);
            x = dx:dx:(n * dx);
            idxOut = [find(sqr' > topThreshold), find(sqr' < bottomThreshold)];
            idxIn =  intersect(find(sqr' < topThreshold), find(sqr' > bottomThreshold)); 
            if n > 1000
                plot(fSQR, x(idxOut), sqr(idxOut), '.', 'color', redColor);
                plot(fSQR, x(idxIn), sqr(idxIn), '.', 'color', blueColor);
            else
                plot(fSQR, x(idxOut), sqr(idxOut), '.', 'color', redColor, 'MarkerSize', 10);
                plot(fSQR, x(idxIn), sqr(idxIn), '.', 'color', blueColor, 'MarkerSize', 10); 
            end
            xlabel('Mean');
            ylabel('SQR');
        end
        
        if plotBoth       
            hBoth = figure;
            fBoth = axes('parent', hBoth);
            hold on;
            yyaxis(fBoth, 'right');
            set(fBoth, 'ycolor', blueColor);
            
            [topThreshold, bottomThreshold] = PlotBeta(data, true);
            idxOut = [find(sqr' > topThreshold), find(sqr' < bottomThreshold)];
            idxIn =  intersect(find(sqr' < topThreshold), find(sqr' > bottomThreshold)); 
            ylim(fBoth, [-2 2]);
            plot([min(data), max(data)], [-2, -2], 'color', 'black');
            ylabel('SQR');
                     
            plot(fBoth, data(idxOut), sqr(idxOut), '.', 'color', redColor);
            plot(fBoth, data(idxIn), sqr(idxIn), '.', 'color', blueColor);                   
            yyaxis(fBoth, 'left');
            plot(fBoth, y, py, '-', 'color', 'black');                     
                        
            [~, idxUnique, ~] = unique(py);
            if length(idxUnique) > 1
                sampleOut = interp1(y(idxUnique), py(idxUnique), data(idxOut));
                plot(fBoth, data(idxOut), sampleOut, '.', 'color', redColor);
            end
            
            set(gca, 'SortMethod', 'depth');                         
            fBoth.YColor = 'k';
            ylabel('PDF');                 
            xlabel('Sample Range');
        end        
    end
end
