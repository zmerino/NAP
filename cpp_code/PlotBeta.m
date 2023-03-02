function [a99, b99] = plotBeta(samples, scaledBySample)
% plotBeta produces a plot of a shaded background, with darker shading
% representing regions of higher probability based on position-dependent
% beta distributions according to sort order statistics.
%
% [a99, b99] = plotBeta(X, SCALED) creates a background plot for length(X)
% positions along the x-axis.  If SCALED is true, the background is scaled
% along the x-axis according to the sample data in X.  if SCALED is false,
% the x-axis is uniform, and dashed contour lines are plotted that represent
% a 99% confidence threshold. These two contour lines are returned as
% output.

    n = length(samples);
    nOriginal = n;
    [a99, ~, a75, ~, b75, ~, b99] = GetTargets(n);   
    skip = 1;
    maxN = 600;
    if n > maxN
        skip = ceil(n / maxN);
        n = ceil(n / skip);
    end
    a = 1:n;
    meanBeta = a ./ (n + 1);          
  
    scale = sqrt(n + 2);
    jInc = (scale) / (n - 1);
    j = -scale:2*jInc:scale;

    yInc = 1 / (n - 1);
    yRange = 0:yInc:1;
    densityMap = zeros(length(yRange), n);
    rangeMap = zeros(length(yRange), n);

    for i = 1:n
        pdf = betapdf(yRange, a(i), n - a(i) + 1);
        rangeR = (yRange - meanBeta(i)) .* scale;
        densityMap(:, i) = pdf;
        rangeMap(:, i) = rangeR;
    end

    c = flipud(gray);
    colormap(c(1:30, :));
    densityMap = densityMap ./ max(densityMap);

    if scaledBySample
        indices = 1:skip:nOriginal; 
        plotRange = samples(indices);
        [X, ~] = meshgrid(plotRange, j);        
        pcolor(X, rangeMap, densityMap);
    else       
        [X, ~] = meshgrid(meanBeta, j);   
        pcolor(X, rangeMap, densityMap);        
        x = (1:nOriginal) / (nOriginal + 1);
        indices = 1:skip:nOriginal;     
        
        hold on;
        g = 128;
        grayScale = [g/255, g/255, g/255];
        plot(x(indices), a99(indices), '--', 'color', grayScale);
        plot(x(indices), b99(indices), '--', 'color', grayScale);

        g = 0;
        grayScale = [g/255, g/255, g/255];
        plot(x(indices), a75(indices), '--', 'color', grayScale);
        plot(x(indices), b75(indices), '--', 'color', grayScale);

        axis([0 1 -2 2])
        plot([0, 1], [-2, -2], 'color', 'black');
        plot([0, 1], [2, 2], 'color', 'black');
    end
    shading interp;
end
