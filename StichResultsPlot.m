function StichResultsPlot(plotQQandSQR,uref,u,msgModelType,Ns,prefix,saveFIG,sqr)
if plotQQandSQR
    figure('Name','SQR and QQ Plots')
    subplot(2,1,1)
    hold off;
    plot(uref,u,'-k','linewidth',1.0);
    xlabel('exact quantile');
    ylabel('empirical quantile');
    title(['QQ-plot:  ',msgModelType,'  N_s = ',num2str(Ns)]);
    xlim([0,1])
    
    subplot(2,1,2)
    hold on;
    % ----------------------------------- create lemon drop oval in gray scale
    smallN = 256;
    smallN2 = 258;
    graymax = 220;
    range = 0:1/(smallN+1):1;
    muLD = range*(smallN + 1) / (smallN + 1);
    lemonDrop = sqrt(muLD.*(1-muLD)) * 3.4;
    sampleCount2 = (smallN + 2):-1:1;
    colorRange = (255-graymax)*sampleCount2/(smallN + 2);
    base = repmat(graymax, smallN + 2, 1);
    col = (base + colorRange') / 255;
    rgb = [col col col];
    count2 = 1;
    
    for ii = ceil(smallN2/2):smallN2-1
        ix = [ii ii+1 smallN2-ii smallN2-ii+1];
        fill(range(ix), lemonDrop(ix), rgb(count2, :),'edgecolor','none')
        fill(range(ix), -lemonDrop(ix), rgb(count2, :),'edgecolor','none')
        count2 = count2 + 2;
    end
    hb1 = plot(muLD,lemonDrop,'k--');
    hb2 = plot(muLD,-lemonDrop,'k--');
    % ------------------------------------------------------------------------
    plot(u,sqr,'-k');
    xlabel('exact quantile');
    ylabel('SQR');
    title(['SQR-plot:  ',msgModelType,'  N_s = ',num2str(Ns)]);
    if saveFIG
        fig5Name = [prefix,'_QQSQR'];
        savefig(fig5Name);
%         print(fig5Name,'-dpng')
    end
    
end
end