classdef distributions
    properties
        x
        Ns
        filename
        pdf_y
        cdf_y
        actual_data
        rndData
        true_quantile
        smooth = 0.01;
        Ncdf = 1000000;
        Npdf = 1000;
        Nhist = 10000;
        precision = 15;
        min_limit = 0;
        max_limit = 10;
        generate_data = false;
        pVector = linspace(0,1,10000);
        dist_name
        distInfo
        randomVSactual = "actual"
        distributionList = ["Beta-a0p5-b1p5";"Beta-a2-b0p5";...
            "Beta-a0p5-b0p5";"Bimodal-Normal";"BirnbaumSaunders";...
            "BirnbaumSaunders-Stable";"Burr";"Exponential";...
            "Extreme-Value";"Gamma";"Generalized-Extreme-Value";...
            "Generalized-Pareto";"HalfNormal";"Normal";...
            "Square-periodic";"Stable";"Stable1";"Stable2";...
            "Stable3";"tLocationScale";"Uniform";"Uniform-Mix";...
            "Weibull";"Chisquare";"InverseGaussian";...
            "Trimodal-Normal"];
    end
    methods
        function [obj] = dist_list(obj)
            debug = false;
            switch obj.dist_name
                case 'Beta-a0p5-b1p5'
                    % Beta1 Case Statement
                    % First shape obj
                    a = 0.5;
                    % Second shape obj
                    b = 1.5;
                    % PDF Curve \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
                    obj.distInfo = makedist('Beta','a',a,'b',b);
                    obj.pdf_y = pdf(obj.distInfo,obj.x);
                    obj.rndData = random(obj.distInfo,1,obj.Ns);
                    obj.actual_data = vertcat(obj.x,obj.pdf_y);
                    % CDF and Q Curve \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
                    obj.cdf_y = cdf(obj.distInfo,obj.x);
                    obj.true_quantile = icdf(obj.distInfo,obj.pVector);
                    %\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
                case 'Beta-a2-b0p5'
                    % Beta2 Case Statement
                    % First shape obj
                    a = 2;
                    % Second shape obj
                    b = 0.5;
                    % PDF Curve \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
                    obj.distInfo = makedist('Beta','a',a,'b',b);
                    obj.pdf_y = pdf(obj.distInfo,obj.x);
                    obj.rndData = random(obj.distInfo,1,obj.Ns);
                    obj.actual_data = vertcat(obj.x,obj.pdf_y);
                    % CDF and Q Curve \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
                    obj.cdf_y = cdf(obj.distInfo,obj.x);
                    obj.true_quantile = icdf(obj.distInfo,obj.pVector);
                    %\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
                case 'Beta-a0p5-b0p5'
                    % Beta3 Case Statement
                    % First shape obj
                    a = 0.5;
                    % Second shape obj
                    b = 0.5;
                    % PDF Curve \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
                    obj.distInfo = makedist('Beta','a',a,'b',b);
                    obj.pdf_y = pdf(obj.distInfo,obj.x);
                    obj.rndData = random(obj.distInfo,1,obj.Ns);
                    obj.actual_data = vertcat(obj.x,obj.pdf_y);
                    % CDF and Q Curve \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
                    obj.cdf_y = cdf(obj.distInfo,obj.x);
                    obj.true_quantile = icdf(obj.distInfo,obj.pVector);
                    %\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
                case 'Bimodal-Normal' %------------------------------------------------- Mixture****
                    % Normal Case Statement
                    % mixture weights
                    p1 = 0.65;
                    p2 = 1 - p1;
                    p = [p1,p2];
                    % Mean
                    Mu1 = 2;
                    Mu2 = 6;
                    % Standard deviation
                    Sigma1 = 0.8;
                    Sigma2 = 0.3;
                    % PDF Curve \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
                    % Distribution 1
                    distributionLabel1 = 'Normal';
                    distInfo1 = makedist(distributionLabel1,...
                        'Mu', Mu1, 'Sigma', Sigma1);
                    pdfCurve1 = pdf(distInfo1,obj.x);
                    % Distribution 2
                    distInfo2 = makedist(distributionLabel1,...
                        'Mu', Mu2, 'Sigma', Sigma2);
                    pdfCurve2 = pdf(distInfo2,obj.x);
                    % Mixture PDF Curve
                    obj.pdf_y = p(1)*pdfCurve1 + p(2)*pdfCurve2;
                    obj.actual_data = vertcat(obj.x,obj.pdf_y);
                    %\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
                    % generate random sample or actual pdf
                    if obj.randomVSactual == "random"
                        % mixture string array flag for mixSampling()
                        mixtureType = "two";
                        % generate n vector for mixture samplings
                        n = misc_functions.mixSampling(obj.Ns,p,mixtureType);
                        % generate random sample
                        rndData1 = random(distInfo1,1,n(1));
                        rndData2 = random(distInfo2,1,n(2));
                        obj.rndData = [rndData1,rndData2];
                    elseif obj.randomVSactual == "actual"
                        obj.actual_data = vertcat(obj.x,obj.pdf_y);
                    end
                    % CREATE DISTRIBUTION OBJECT --------------------------
                    % mixture string array flag for misc_functions.mixSampling()
                    mixtureType = "two";
                    % generate m vector for mixture samplings
                    m = misc_functions.mixSampling(obj.Ncdf,p,mixtureType);
                    % generate random sample to create distribution
                    % object
                    actData1 = random(distInfo1,1,m(1));
                    actData2 = random(distInfo2,1,m(2));
                    actual_data = [actData1,actData2];
                    % generate numerical cdf: f=cdf, s=x-coordinates
                    [f,s] = ecdf(actual_data);
                    f = f(1:2:end);
                    s = s(1:2:end);
                    % generate distribution object
                    obj.distInfo = ...
                        makedist('PiecewiseLinear','x',s,'Fx',f);
                    obj.cdf_y = cdf(obj.distInfo,obj.x);
                    % create cdf/pdf from distribution object. 
                    % for debugging and vizualization.
                    if debug
                        if obj.randomVSactual == "actual"
                            xMix = linspace(obj.min_limit,...
                                obj.max_limit,obj.Npdf);
                            CDF = cdf(obj.distInfo,xMix);
                            % numerically differentiate CDF
                            PDF = zeros(1,size(CDF(1:end-1),2));
                            for i = 1:size(CDF,2)-1
                                dx1 = (xMix(i+1)-xMix(i));
                                PDF(i) =  (CDF(i+1) - CDF(i))/dx1;
                            end
                            % smooth pdf obj.actual_data
                            smoo1 = smooth(xMix(1:end-1),PDF,obj.smooth);
                            % plot cdf,pdf,smoothed-pdf
                            figure('Name',['Debug: ',...
                                char(obj.dist_name)])
                            subplot(2,1,1)
                            hold on
                            plot(xMix,CDF,'-m')
                            plot(xMix(1:end-1),PDF,'-r')
                            plot(xMix(1:end-1),smoo1,'-b')
                            ylabel('$f(x) or F(x)$','Interpreter','latex')
                            xlabel('x','Interpreter','latex')
                            legend('cdf','pdf','smoothed-pdf')
                            % plot histogram for random sample
                            subplot(2,1,2)
                            histogram(random(obj.distInfo,obj.Nhist,1),...
                                'Normalization','probability')
                            xlim([obj.min_limit,obj.max_limit])
                        end
                    end
                    %------------------------------------------------------
                case 'Binomial'
                    % Binomial Case Statement
                    % Number of trials
                    n = 2000;
                    % Porbability of success for each trial
                    p = 0.2;
                    % PDF Curve \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
                    obj.distInfo = ...
                        makedist(obj.dist_name,'n',n,'p',p);
                    obj.pdf_y = binopdf(obj.x,n,p);
                    obj.rndData = random(obj.distInfo,1,obj.Ns);
                    obj.actual_data = vertcat(obj.x,obj.pdf_y);
                    % CDF and Q Curve \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
                    obj.cdf_y = cdf(obj.distInfo,obj.x);
                    obj.true_quantile = icdf(obj.distInfo,obj.pVector);
                    %\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
                case 'BirnbaumSaunders'
                    % BirnbaumSaunders Case Statement
                    % Scale parameter
                    Beta = 1.5;
                    % Shape parameter
                    Gamma = 0.5;
                    % PDF Curve \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
                    obj.distInfo = makedist(obj.dist_name,...
                        'Beta',Beta,'Gamma',Gamma);
                    obj.pdf_y = pdf(obj.distInfo,obj.x);
                    obj.rndData = random(obj.distInfo,1,obj.Ns);
                    obj.actual_data = vertcat(obj.x,obj.pdf_y);
                    % CDF and Q Curve \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
                    obj.cdf_y = cdf(obj.distInfo,obj.x);
                    obj.true_quantile = icdf(obj.distInfo,obj.pVector);
                    %\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
                case 'BirnbaumSaunders-Stable' %------------------------------------------------- Mixture****
                    % BirnbaumSaunders Case Statement
                    % mixture weights
                    p1 = 0.35;
                    p2 = 1 - p1;
                    p = [p1,p2];
                    % BirnbaumSaunders distribution -----------------------
                    % Scale parameter
                    Beta = 1.5;
                    % Shape parameter
                    Gamma = 0.5;
                    % Stable distribution ---------------------------------
                    % First shape parameter
                    Alpha1 = 0.5;
                    % Second shape parameter: -1 <= Beta <= 1
                    Beta1 = 0.05;
                    % Scale parameter
                    Gam1 = 1;
                    % Location parameter
                    Delta1 = 7;
                    % PDF Curve \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
                    % BirnbaumSaunders distribution
                    distributionLabel1 = 'BirnbaumSaunders';
                    distInfo1 = makedist(distributionLabel1,...
                        'Beta',Beta,'Gamma',Gamma);
                    pdfCurve1 = pdf(distInfo1,obj.x);
                    % Stable distribution
                    distributionLabel2 = 'Stable';
                    distInfo2 = makedist(distributionLabel2,...
                        'Alpha', Alpha1,'Beta', Beta1,...
                        'Gam', Gam1, 'Delta', Delta1);
                    pdfCurve2 = pdf(distInfo2, obj.x);
                    % Mixture PDF Curve
                    obj.pdf_y = p(1)*pdfCurve1 + p(2)*pdfCurve2;
                    obj.actual_data = vertcat(obj.x,obj.pdf_y);
                    %\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
                    % generate random sample or actual pdf
                    if obj.randomVSactual == "random"
                        % mixture string array flag for misc_functions.mixSampling()
                        mixtureType = "two";
                        % generate n vector for mixture samplings
                        n = misc_functions.mixSampling(obj.Ns,p,mixtureType);
                        % generate random sample
                        rndData1 = random(distInfo1,1,n(1));
                        rndData2 = random(distInfo2,1,n(2));
                        obj.rndData = [rndData1,rndData2];
                    elseif obj.randomVSactual == "actual"
                        obj.actual_data = vertcat(obj.x,obj.pdf_y);
                    end
                    % CREATE DISTRIBUTION OBJECT --------------------------
                    % mixture string array flag for misc_functions.mixSampling()
                    mixtureType = "two";
                    % generate m vector for mixture samplings
                    m = misc_functions.mixSampling(obj.Ncdf,p,mixtureType);
                    % generate random sample to create distribution
                    % object
                    actData1 = random(distInfo1,1,m(1));
                    actData2 = random(distInfo2,1,m(2));
                    obj.actual_data = [actData1,actData2];
                    % generate numerical cdf: f=cdf, s=x-coordinates
                    [f,s] = ecdf(obj.actual_data);
                    f = f(1:2:end);
                    s = s(1:2:end);
                    % generate distribution object
                    obj.distInfo = ...
                        makedist('PiecewiseLinear','x',s,'Fx',f);
                    obj.cdf_y = cdf(obj.distInfo,obj.x);
                    % create cdf/pdf from distribution object. 
                    % for debugging and vizualization.
                    if debug
                        if obj.randomVSactual == "actual"
                            xMix = linspace(obj.min_limit,...
                                obj.max_limit,obj.Npdf);
                            CDF = cdf(obj.distInfo,xMix);
                            % numerically differentiate
                            PDF = zeros(1,size(CDF(1:end-1),2));
                            for i = 1:size(CDF,2)-1
                                dx1 = (xMix(i+1)-xMix(i));
                                PDF(i) =  (CDF(i+1) - CDF(i))/dx1;
                            end
                            % smooth pdf obj.actual_data
                            smoo1 = smooth(xMix(1:end-1),PDF,obj.smooth);
                            % plot cdf,pdf,smoothed-pdf
                            figure('Name',['Debug: ',...
                                char(obj.dist_name)])
                            subplot(2,1,1)
                            hold on
                            plot(xMix,CDF,'-m')
                            plot(xMix(1:end-1),PDF,'-r')
                            plot(xMix(1:end-1),smoo1,'-b')
                            ylabel('$f(x) or F(x)$','Interpreter','latex')
                            xlabel('x','Interpreter','latex')
                            legend('cdf','pdf','smoothed-pdf')
                            % plot histogram for random sample
                            subplot(2,1,2)
                            histogram(random(obj.distInfo,obj.Nhist,1),...
                                'Normalization','probability')
                            xlim([obj.min_limit,obj.max_limit])
                        end
                    end
                    %------------------------------------------------------
                case 'Burr'
                    % Burr Case Statement
                    % Scale parameter
                    Alpha = 1;
                    % Shape parameter one
                    c = 2;
                    % Shape parameter two
                    k = 2;
                    % PDF Curve \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
                    obj.distInfo = makedist(obj.dist_name,...
                        'Alpha',Alpha,'c',c,'k',k);
                    obj.pdf_y = pdf(obj.distInfo,obj.x);
                    obj.rndData = random(obj.distInfo,1,obj.Ns);
                    obj.actual_data = vertcat(obj.x,obj.pdf_y);
                    % CDF and Q Curve \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
                    obj.cdf_y = cdf(obj.distInfo,obj.x);
                    obj.true_quantile = icdf(obj.distInfo,obj.pVector);
                    %\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
                case 'Chisquare'
                    % Chisquare Case Statement
                    % Degrees of freedom
                    Nu = 4;
                    % PDF Curve \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
                    obj.pdf_y = chi2pdf(obj.x,Nu);
                    obj.rndData = random(obj.distInfo,1,obj.Ns);
                    obj.actual_data = vertcat(obj.x,obj.pdf_y);
                    % CDF and Q Curve \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
                    obj.cdf_y = cdf(obj.distInfo,obj.x);
                    obj.true_quantile = icdf(obj.distInfo,obj.pVector);
                    %\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
                case 'Exponential'
                    % Exponential Case Statement
                    % Mean
                    Mu = 1;
                    % PDF Curve \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
                    obj.distInfo = makedist(obj.dist_name,'Mu',Mu);
                    obj.pdf_y = pdf(obj.distInfo,obj.x);
                    obj.rndData = random(obj.distInfo,1,obj.Ns);
                    obj.actual_data = vertcat(obj.x,obj.pdf_y);
                    % CDF and Q Curve \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
                    obj.cdf_y = cdf(obj.distInfo,obj.x);
                    obj.true_quantile = icdf(obj.distInfo,obj.pVector);
                    %\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
                case 'Extreme-Value'
                    % Extreme Value Case Statement
                    % Location parameter
                    Mu = 1;
                    % Scale parameter
                    Sigma = 2;
                    % PDF Curve \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
                    distributionLabel = 'Extreme Value';
                    obj.distInfo = makedist(distributionLabel,...
                        'Mu',Mu, 'Sigma', Sigma);
                    obj.pdf_y = pdf(obj.distInfo,obj.x);
                    obj.rndData = random(obj.distInfo,1,obj.Ns);
                    obj.actual_data = vertcat(obj.x,obj.pdf_y);
                    % CDF and Q Curve \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
                    obj.cdf_y = cdf(obj.distInfo,obj.x);
                    obj.true_quantile = icdf(obj.distInfo,obj.pVector);
                    %\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
                case 'Gamma'
                    % Gamma Case Statement
                    % Shape parameter
                    a = 2;
                    % Scale parameter
                    b = 2;
                    % PDF Curve \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
                    obj.distInfo = ...
                        makedist(obj.dist_name,'a',a,'b',b);
                    obj.pdf_y = pdf(obj.distInfo,obj.x);
                    obj.rndData = random(obj.distInfo,1,obj.Ns);
                    obj.actual_data = vertcat(obj.x,obj.pdf_y);
                    % CDF and Q Curve \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
                    obj.cdf_y = cdf(obj.distInfo,obj.x);
                    obj.true_quantile = icdf(obj.distInfo,obj.pVector);
                    %\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
                case 'Generalized-Extreme-Value' %------------------------------------------------- Mixture****
                    % Generalized Extreme Value Value Statement
                    % Shape parameter
                    k = 1;
                    % Scale parameter
                    Sigma = 2;
                    % Location parameter
                    Mu = 2;
                    % PDF Curve \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
                    distributionLabel = 'Generalized Extreme Value';
                    obj.distInfo = makedist(distributionLabel,...
                        'k',k, 'Sigma', Sigma,'Mu',Mu);
                    obj.pdf_y = pdf(obj.distInfo,obj.x);
                    obj.rndData = random(obj.distInfo,1,obj.Ns);
                    obj.actual_data = vertcat(obj.x,obj.pdf_y);
                    % CDF and Q Curve \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
                    obj.cdf_y = cdf(obj.distInfo,obj.x);
                    obj.true_quantile = icdf(obj.distInfo,obj.pVector);
                    %\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
                case 'Generalized-Pareto' %------------------------------------------------- Mixture****
                    % Generalized Pareto Value Value Case Statement
                    % Tail MemTracker (shape) parameter
                    k = 2;
                    % Scale parameter
                    Sigma = 1;
                    % Threshold (location) parameter
                    theta = 0;
                    % PDF Curve \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
                    distributionLabel = 'Generalized Pareto';
                    obj.distInfo = makedist(distributionLabel,...
                        'k',k, 'Sigma', Sigma,'Theta', theta);
                    obj.pdf_y = pdf(obj.distInfo,obj.x);
                    % Mu parameter is not recognized
                    obj.rndData = random(obj.distInfo,1,obj.Ns);
                    obj.actual_data = vertcat(obj.x,obj.pdf_y);
                    % CDF and Q Curve \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
                    obj.cdf_y = cdf(obj.distInfo,obj.x);
                    obj.true_quantile = icdf(obj.distInfo,obj.pVector);
                    %\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
                case 'HalfNormal'
                    % Half Normal Value Case Statement
                    % Location parameter
                    Mu = 0;
                    % Scale parameter
                    Sigma = 1;
                    % PDF Curve \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
                    obj.distInfo = makedist(obj.dist_name,...
                        'Mu', Mu, 'Sigma', Sigma);
                    obj.pdf_y = pdf(obj.dist_name,obj.x);
                    obj.rndData = random(obj.distInfo,1,obj.Ns);
                    obj.actual_data = vertcat(obj.x,obj.pdf_y);
                    % CDF and Q Curve \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
                    obj.cdf_y = cdf(obj.distInfo,obj.x);
                    obj.true_quantile = icdf(obj.distInfo,obj.pVector);
                    %\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
                case 'InverseGaussian'
                    % Inverse Gaussian Case Statement
                    % Scale parameter
                    Mu = 1;
                    % Shape parameter
                    Lambda = 1;
                    % PDF Curve \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
                    obj.distInfo = makedist(obj.dist_name,...
                        'mu', Mu, 'lambda', Lambda);
                    obj.pdf_y = pdf(obj.distInfo,obj.x);
                    obj.rndData = random(obj.distInfo,1,obj.Ns);
                    obj.actual_data = vertcat(obj.x,obj.pdf_y);
                    % CDF and Q Curve \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
                    obj.cdf_y = cdf(obj.distInfo,obj.x);
                    obj.true_quantile = icdf(obj.distInfo,obj.pVector);
                    %\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
                case 'Normal'
                    % Normal Case Statement
                    % Mean
                    Mu = 5;
                    % Standard deviation
                    Sigma = 1;
                    % PDF Curve \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
                    obj.distInfo = makedist(obj.dist_name,...
                        'Mu', Mu, 'Sigma', Sigma);
                    obj.pdf_y = pdf(obj.distInfo,obj.x);
                    obj.rndData = random(obj.distInfo,1,obj.Ns);
                    obj.actual_data = vertcat(obj.x,obj.pdf_y);
                    % CDF and Q Curve \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
                    obj.cdf_y = cdf(obj.distInfo,obj.x);
                    obj.true_quantile = icdf(obj.distInfo,obj.pVector);
                    %\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
                case 'Normal-Contaminated' %------------------------------------------------- Mixture****
                    % Normal Case Statement
                    % mixture weights
                    p1 = 0.5;
                    p2 = 1 - p1;
                    p = [p1,p2];
                    % Mean
                    Mu1 = 5;
                    Mu2 = 5;
                    % Standard deviation
                    Sigma1 = 2;
                    Sigma2 = 0.25;
                    % PDF Curve \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
                    % Distribution 1
                    distributionLabel1 = 'Normal';
                    distInfo1 = makedist(distributionLabel1,...
                        'Mu', Mu1, 'Sigma', Sigma1);
                    pdfCurve1 = pdf(distInfo1,obj.x);
                    % Distribution 2
                    distInfo2 = makedist(distributionLabel1,...
                        'Mu', Mu2, 'Sigma', Sigma2);
                    pdfCurve2 = pdf(distInfo2,obj.x);
                    % Mixture PDF Curve
                    obj.pdf_y = p(1)*pdfCurve1 + p(2)*pdfCurve2;
                    obj.actual_data = vertcat(obj.x,obj.pdf_y);
                    %\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
                    % generate random sample or actual pdf
                    if obj.randomVSactual == "random"
                        % mixture string array flag for misc_functions.mixSampling()
                        mixtureType = "two";
                        % generate n vector for mixture samplings
                        n = misc_functions.mixSampling(obj.Ns,p,mixtureType);
                        % generate random sample
                        rndData1 = random(distInfo1,1,n(1));
                        rndData2 = random(distInfo2,1,n(2));
                        obj.rndData = [rndData1,rndData2];
                    elseif obj.randomVSactual == "actual"
                        obj.actual_data = vertcat(obj.x,obj.pdf_y);
                    end
                    % CREATE DISTRIBUTION OBJECT --------------------------
                    % mixture string array flag for misc_functions.mixSampling()
                    mixtureType = "two";
                    % generate m vector for mixture samplings
                    m = misc_functions.mixSampling(obj.Ncdf,p,mixtureType);
                    % generate random sample to create distribution
                    % object
                    actData1 = random(distInfo1,1,m(1));
                    actData2 = random(distInfo2,1,m(2));
                    obj.actual_data = [actData1,actData2];
                    % generate numerical cdf: f=cdf, s=x-coordinates
                    [f,s] = ecdf(obj.actual_data);
                    f = f(1:2:end);
                    s = s(1:2:end);
                    % generate distribution object
                    obj.distInfo = ...
                        makedist('PiecewiseLinear','x',s,'Fx',f);
                    obj.cdf_y = cdf(obj.distInfo,obj.x);
                    % create cdf/pdf from distribution object. 
                    % for debugging and vizualization.
                    if debug
                        if obj.randomVSactual == "actual"
                            xMix = linspace(obj.min_limit,...
                                obj.max_limit,obj.Npdf);
                            CDF = cdf(obj.distInfo,xMix);
                            % numerically differentiate
                            PDF = zeros(1,size(CDF(1:end-1),2));
                            for i = 1:size(CDF,2)-1
                                dx1 = (xMix(i+1)-xMix(i));
                                PDF(i) =  (CDF(i+1) - CDF(i))/dx1;
                            end
                            % smooth pdf obj.actual_data
                            smoo1 = smooth(xMix(1:end-1),PDF,obj.smooth);
                            % plot cdf,pdf,smoothed-pdf
                            figure('Name',['Debug: ',...
                                char(obj.dist_name)])
                            subplot(2,1,1)
                            hold on
                            plot(xMix,CDF,'-m')
                            plot(xMix(1:end-1),PDF,'-r')
                            plot(xMix(1:end-1),smoo1,'-b')
                            ylabel('$f(x) or F(x)$','Interpreter','latex')
                            xlabel('x','Interpreter','latex')
                            legend('cdf','pdf','smoothed-pdf')
                            % plot histogram for random sample
                            subplot(2,1,2)
                            histogram(random(obj.distInfo,obj.Nhist,1),...
                                'Normalization','probability')
                            xlim([obj.min_limit,obj.max_limit])
                        end
                    end
                    %------------------------------------------------------
                case 'Square-periodic' %------------------------------------------------- Mixture****
                    % Uniform Case Statement
                    % mixture weights
                    p1 = 1/6;
                    p2 = 1/6;
                    p3 = 1/6;
                    p4 = 1/6;
                    p5 = 1/6;
                    p6 = 1 - p1 - p2 - p3 - p4- p5 ;
                    p = [p1,p2,p3,p4,p5,p6];
                    % Lower bound
                    Lower1 = 1;
                    Lower2 = 2.5;
                    Lower3 = 4;
                    Lower4 = 5.5;
                    Lower5 = 7;
                    Lower6 = 8.5;
                    % Upper Bound
                    Upper1 = 2;
                    Upper2 = 3.5;
                    Upper3 = 5;
                    Upper4 = 6.5;
                    Upper5 = 8;
                    Upper6 = 9.5;
                    % PDF Curve \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
                    distributionLabel1 = 'Uniform';
                    % Distribution 1
                    distInfo1 = makedist(distributionLabel1,...
                        'Lower', Lower1, 'Upper', Upper1);
                    pdfCurve1 = pdf(distInfo1,obj.x);
                    % Distribution 2
                    distInfo2 = makedist(distributionLabel1,...
                        'Lower', Lower2, 'Upper', Upper2);
                    pdfCurve2 = pdf(distInfo2,obj.x);
                    % Distribution 3
                    distInfo3 = makedist(distributionLabel1,...
                        'Lower', Lower3, 'Upper', Upper3);
                    pdfCurve3 = pdf(distInfo3,obj.x);
                    % Distribution 4
                    distInfo4 = makedist(distributionLabel1,...
                        'Lower', Lower4, 'Upper', Upper4);
                    pdfCurve4 = pdf(distInfo4,obj.x);
                    % Distribution 5
                    distInfo5 = makedist(distributionLabel1,...
                        'Lower', Lower5, 'Upper', Upper5);
                    pdfCurve5 = pdf(distInfo5,obj.x);
                    % Distribution 6
                    distInfo6 = makedist(distributionLabel1,...
                        'Lower', Lower6, 'Upper', Upper6);
                    pdfCurve6 = pdf(distInfo6,obj.x);
                    % Mixture PDF Curve
                    obj.pdf_y = p(1)*pdfCurve1 + p(2)*pdfCurve2 +...
                        p(3)*pdfCurve3 + p(4)*pdfCurve4 +...
                        p(5)*pdfCurve5 + p(6)*pdfCurve6;
                    %\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
                    % generate random sample or actual pdf
                    if obj.randomVSactual == "random"
                        % mixture string array flag for misc_functions.mixSampling()
                        mixtureType = "six";
                        % generate n vector for mixture samplings
                        n = misc_functions.mixSampling(obj.Ns,p,mixtureType);
                        % generate random sample
                        rndData1 = random(distInfo1,1,n(1));
                        rndData2 = random(distInfo2,1,n(2));
                        rndData3 = random(distInfo3,1,n(3));
                        rndData4 = random(distInfo4,1,n(4));
                        rndData5 = random(distInfo5,1,n(5));
                        rndData6 = random(distInfo6,1,n(6));
                        obj.rndData = [rndData1,rndData2,rndData3,...
                            rndData4,rndData5,rndData6];
                    elseif obj.randomVSactual == "actual"
                        obj.actual_data = vertcat(obj.x,obj.pdf_y);
                    end
                    % CREATE DISTRIBUTION OBJECT --------------------------
                    % mixture string array flag for misc_functions.mixSampling()
                    mixtureType = "six";
                    % generate m vector for mixture samplings
                    m = misc_functions.mixSampling(obj.Ncdf,p,mixtureType);
                    % generate random sample to create distribution
                    % object
                    actData1 = random(distInfo1,1,m(1));
                    actData2 = random(distInfo2,1,m(2));
                    actData3 = random(distInfo3,1,m(3));
                    actData4 = random(distInfo4,1,m(4));
                    actData5 = random(distInfo5,1,m(5));
                    actData6 = random(distInfo6,1,m(6));
                    obj.actual_data = [actData1,actData2,actData3,...
                        actData4,actData5,actData6];
                    % generate numerical cdf: f=cdf, s=x-coordinates
                    [f,s] = ecdf(obj.actual_data);
                    f = f(1:2:end);
                    s = s(1:2:end);
                    % generate distribution object
                    obj.distInfo = ...
                        makedist('PiecewiseLinear','x',s,'Fx',f);
                    obj.cdf_y = cdf(obj.distInfo,obj.x);
                    % create cdf/pdf from distribution object. 
                    % for debugging and vizualization.
                    if debug
                        if obj.randomVSactual == "actual"
                            xMix = linspace(obj.min_limit,...
                                obj.max_limit,obj.Npdf);
                            CDF = cdf(obj.distInfo,xMix);
                            % numerically differentiate
                            PDF = zeros(1,size(CDF(1:end-1),2));
                            for i = 1:size(CDF,2)-1
                                dx1 = (xMix(i+1)-xMix(i));
                                PDF(i) =  (CDF(i+1) - CDF(i))/dx1;
                            end
                            % smooth pdf obj.actual_data
                            smoo1 = smooth(xMix(1:end-1),PDF,obj.smooth);
                            % plot cdf,pdf,smoothed-pdf
                            figure('Name',['Debug: ',...
                                char(obj.dist_name)])
                            subplot(2,1,1)
                            hold on
                            plot(xMix,CDF,'-m')
                            plot(xMix(1:end-1),PDF,'-r')
                            plot(xMix(1:end-1),smoo1,'-b')
                            ylabel('$f(x) or F(x)$','Interpreter','latex')
                            xlabel('x','Interpreter','latex')
                            legend('cdf','pdf','smoothed-pdf')
                            % plot histogram for random sample
                            subplot(2,1,2)
                            histogram(random(obj.distInfo,obj.Nhist,1),...
                                'Normalization','probability')
                            xlim([obj.min_limit,obj.max_limit])
                        end
                    end
                    %}
                    %------------------------------------------------------
                case 'Stable'
                    % Stable Case Statement
                    % First shape parameter
                    Alpha = 0.5;
                    %{
                    Alpha = 0.4;
                    Alpha = 0.35
                    Alpha = 0.5  %<---- use to start
                    Alpha = 0.2; %<---- very hard to estimate
                    %}
                    % Second shape parameter: -1 <= Beta <= 1
                    Beta = 0.05;
                    %{
                    Beta = 0.9;
                    Beta = 1;
                    Beta = .05;
                    Beta = .05;
                    %}
                    % Scale parameter
                    Gam = 1;
                    % Location parameter
                    Delta = 4;
                    % PDF Curve \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
                    obj.distInfo = makedist(obj.dist_name,...
                        'Alpha', Alpha,'Beta', Beta,...
                        'Gam', Gam, 'Delta', Delta);
                    obj.pdf_y = pdf(obj.distInfo, obj.x);
                    obj.rndData = random(obj.distInfo,1,obj.Ns);
                    obj.actual_data = vertcat(obj.x,obj.pdf_y);
                    % CDF and Q Curve \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
                    obj.cdf_y = cdf(obj.distInfo,obj.x);
                    obj.true_quantile = icdf(obj.distInfo,obj.pVector);
                    %\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
                case 'Stable1'
                    % Stable Case Statement
                    % First shape parameter
                    Alpha = 0.2;
                    %{
                    Alpha = 0.4;
                    Alpha = 0.35
                    Alpha = 0.5  %<---- use to start
                    Alpha = 0.2; %<---- very hard to estimate
                    %}
                    % Second shape parameter: -1 <= Beta <= 1
                    Beta = 0.05;
                    %{
                    Beta = 0.9;
                    Beta = 1;
                    Beta = .05;
                    Beta = .05;
                    %}
                    % Scale parameter
                    Gam = 1;
                    % Location parameter
                    Delta = 4;
                    % PDF Curve \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
                    distributionLabel = 'Stable';
                    obj.distInfo = makedist(distributionLabel,...
                        'Alpha', Alpha,'Beta', Beta,...
                        'Gam', Gam, 'Delta', Delta);
                    obj.pdf_y = pdf(obj.distInfo, obj.x);
                    obj.rndData = random(obj.distInfo,1,obj.Ns);
                    obj.actual_data = vertcat(obj.x,obj.pdf_y);
                    % CDF and Q Curve \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
                    obj.cdf_y = cdf(obj.distInfo,obj.x);
                    obj.true_quantile = icdf(obj.distInfo,obj.pVector);
                    %\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
                case 'Stable2' %------------------------------------------------- Mixture****
                    % mixture model for 2 stable distributions
                    % mixture weights
                    p1 = 0.25;
                    p2 = 1 - p1;
                    p = [p1,p2];
                    % First shape parameter
                    Alpha1 = 0.5;
                    Alpha2 = 0.5;
                    % Second shape parameter: -1 <= Beta <= 1
                    Beta1 = 0.05;
                    Beta2 = 0.05;
                    % Scale parameter
                    Gam1 = 1;
                    Gam2 = 1;
                    % Location parameter
                    Delta1 = 2;
                    Delta2 = 5;
                    % PDF Curve \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
                    distributionLabel = 'Stable';
                    % stable 1
                    distInfo1 = makedist(distributionLabel,...
                        'Alpha', Alpha1,'Beta', Beta1,...
                        'Gam', Gam1, 'Delta', Delta1);
                    pdfCurve1 = pdf(distInfo1, obj.x);
                    % stable 2
                    distInfo2 = makedist(distributionLabel,...
                        'Alpha', Alpha2,'Beta', Beta2,...
                        'Gam', Gam2, 'Delta', Delta2);
                    pdfCurve2 = pdf(distInfo2, obj.x);
                    % Mixture PDF Curve
                    obj.pdf_y = p(1)*pdfCurve1 + p(2)*pdfCurve2;
                    %\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
                    % generate random sample or actual pdf
                    if obj.randomVSactual == "random"
                        % mixture string array flag for misc_functions.mixSampling()
                        mixtureType = "two";
                        % generate n vector for mixture samplings
                        n = misc_functions.mixSampling(obj.Ns,p,mixtureType);
                        % generate random sample
                        rndData1 = random(distInfo1,1,n(1));
                        rndData2 = random(distInfo2,1,n(2));
                        obj.rndData = [rndData1,rndData2];
                    elseif obj.randomVSactual == "actual"
                        obj.actual_data = vertcat(obj.x,obj.pdf_y);
                    end
                    % CREATE DISTRIBUTION OBJECT --------------------------
                    % mixture string array flag for misc_functions.mixSampling()
                    mixtureType = "two";
                    % generate m vector for mixture samplings
                    m = misc_functions.mixSampling(obj.Ncdf,p,mixtureType);
                    % generate random sample to create distribution
                    % object
                    actData1 = random(distInfo1,1,m(1));
                    actData2 = random(distInfo2,1,m(2));
                    obj.actual_data = [actData1,actData2];
                    % generate numerical cdf: f=cdf, s=x-coordinates
                    [f,s] = ecdf(obj.actual_data);
                    f = f(1:2:end);
                    s = s(1:2:end);
                    % generate distribution object
                    obj.distInfo = ...
                        makedist('PiecewiseLinear','x',s,'Fx',f);
                    obj.cdf_y = cdf(obj.distInfo,obj.x);
                    % create cdf/pdf from distribution object. 
                    % for debugging and vizualization.
                    if debug
                        if obj.randomVSactual == "actual"
                            xMix = linspace(obj.min_limit,...
                                obj.max_limit,obj.Npdf);
                            CDF = cdf(obj.distInfo,xMix);
                            % numerically differentiate
                            PDF = zeros(1,size(CDF(1:end-1),2));
                            for i = 1:size(CDF,2)-1
                                dx1 = (xMix(i+1)-xMix(i));
                                PDF(i) =  (CDF(i+1) - CDF(i))/dx1;
                            end
                            % smooth pdf obj.actual_data
                            smoo1 = smooth(xMix(1:end-1),PDF,obj.smooth);
                            % plot cdf,pdf,smoothed-pdf
                            figure('Name',['Debug: ',...
                                char(obj.dist_name)])
                            subplot(2,1,1)
                            hold on
                            plot(xMix,CDF,'-m')
                            plot(xMix(1:end-1),PDF,'-r')
                            plot(xMix(1:end-1),smoo1,'-b')
                            ylabel('$f(x) or F(x)$','Interpreter','latex')
                            xlabel('x','Interpreter','latex')
                            legend('cdf','pdf','smoothed-pdf')
                            % plot histogram for random sample
                            subplot(2,1,2)
                            histogram(random(obj.distInfo,obj.Nhist,1),...
                                'Normalization','probability')
                            xlim([obj.min_limit,obj.max_limit])
                        end
                    end
                    %------------------------------------------------------
                case 'Stable3' %------------------------------------------------- Mixture****
                    % mixture model for 3 stable distributions
                    % mixture weights
                    p1 = 0.25;
                    p2 = 0.5;
                    p3 = 1 - p1 - p2;
                    p = [p1,p2,p3];
                    % Stable distributions --------------------------------
                    % First shape parameter
                    Alpha1 = 0.5;
                    Alpha2 = 0.5;
                    Alpha3 = 0.5;
                    % Second shape parameter: -1 <= Beta <= 1
                    Beta1 = 0.05;
                    Beta2 = 0.05;
                    Beta3 = 0.05;
                    % Scale parameter
                    Gam1 = 1;
                    Gam2 = 1;
                    Gam3 = 1;
                    % Location parameter
                    Delta1 = 2;
                    Delta2 = 5;
                    Delta3 = 8;
                    % PDF Curve \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
                    distributionLabel = 'Stable';
                    % stable 1
                    distInfo1 = makedist(distributionLabel,...
                        'Alpha', Alpha1,'Beta', Beta1,...
                        'Gam', Gam1, 'Delta', Delta1);
                    pdfCurve1 = pdf(distInfo1,obj.x);
                    % stable 2
                    distInfo2 = makedist(distributionLabel,...
                        'Alpha', Alpha2,'Beta', Beta2,...
                        'Gam', Gam2, 'Delta', Delta2);
                    pdfCurve2 = pdf(distInfo2,obj.x);
                    % stable 3
                    distInfo3 = makedist(distributionLabel,...
                        'Alpha', Alpha3,'Beta', Beta3,...
                        'Gam', Gam3, 'Delta', Delta3);
                    pdfCurve3 = pdf(distInfo3,obj.x);
                    % Mixture PDF Curve
                    obj.pdf_y = p(1)*pdfCurve1 +...
                        p(2)*pdfCurve2 + p(3)*pdfCurve3;
                    %\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
                    % generate random sample or actual pdf
                    if obj.randomVSactual == "random"
                        % mixture string array flag for misc_functions.mixSampling()
                        mixtureType = "three";
                        % generate n vector for mixture samplings
                        n = misc_functions.mixSampling(obj.Ns,p,mixtureType);
                        % generate random sample
                        rndData1 = random(distInfo1,1,n(1));
                        rndData2 = random(distInfo2,1,n(2));
                        rndData3 = random(distInfo3,1,n(3));
                        obj.rndData = [rndData1,rndData2,rndData3];
                    elseif obj.randomVSactual == "actual"
                        obj.actual_data = vertcat(obj.x,obj.pdf_y);
                    end
                    % CREATE DISTRIBUTION OBJECT --------------------------
                    % mixture string array flag for misc_functions.mixSampling()
                    mixtureType = "three";
                    % generate m vector for mixture samplings
                    m = misc_functions.mixSampling(obj.Ncdf,p,mixtureType);
                    % generate random sample to create distribution
                    % object
                    actData1 = random(distInfo1,1,m(1));
                    actData2 = random(distInfo2,1,m(2));
                    actData3 = random(distInfo3,1,m(3));
                    obj.actual_data = [actData1,actData2,actData3];
                    % generate numerical cdf: f=cdf, s=x-coordinates
                    [f,s] = ecdf(obj.actual_data);
                    f = f(1:2:end);
                    s = s(1:2:end);
                    % generate distribution object
                    obj.distInfo = ...
                        makedist('PiecewiseLinear','x',s,'Fx',f);
                    obj.cdf_y = cdf(obj.distInfo,obj.x);
                    % create cdf/pdf from distribution object. 
                    % for debugging and vizualization.
                    if debug
                        if obj.randomVSactual == "actual"
                            xMix = linspace(obj.min_limit,...
                                obj.max_limit,obj.Npdf);
                            CDF = cdf(obj.distInfo,xMix);
                            % numerically differentiate
                            PDF = zeros(1,size(CDF(1:end-1),2));
                            for i = 1:size(CDF,2)-1
                                dx1 = (xMix(i+1)-xMix(i));
                                PDF(i) =  (CDF(i+1) - CDF(i))/dx1;
                            end
                            % smooth pdf obj.actual_data
                            smoo1 = smooth(xMix(1:end-1),PDF,obj.smooth);
                            % plot cdf,pdf,smoothed-pdf
                            figure('Name',['Debug: ',...
                                char(obj.dist_name)])
                            subplot(2,1,1)
                            hold on
                            plot(xMix,CDF,'-m')
                            plot(xMix(1:end-1),PDF,'-r')
                            plot(xMix(1:end-1),smoo1,'-b')
                            ylabel('$f(x) or F(x)$','Interpreter','latex')
                            xlabel('x','Interpreter','latex')
                            legend('cdf','pdf','smoothed-pdf')
                            % plot histogram for random sample
                            subplot(2,1,2)
                            histogram(random(obj.distInfo,obj.Nhist,1),...
                                'Normalization','probability')
                            xlim([obj.min_limit,obj.max_limit])
                        end
                    end
                    %------------------------------------------------------
                case 'Trimodal-Normal' %------------------------------------------------- Mixture****
                    % Normal Case Statement
                    % mixture weights
                    p1 = 0.33;
                    p2 = 0.33;
                    p3 = 1 - p1 - p2;
                    p = [p1,p2,p3];
                    % Mean
                    Mu1 = 4;
                    Mu2 = 5;
                    Mu3 = 6;
                    % Standard deviation
                    Sigma1 = 0.5;
                    Sigma2 = 0.25;
                    Sigma3 = 0.5;
                    % PDF Curve \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
                    % Distribution 1
                    distributionLabel1 = 'Normal';
                    distInfo1 = makedist(distributionLabel1,...
                        'Mu', Mu1, 'Sigma', Sigma1);
                    pdfCurve1 = pdf(distInfo1,obj.x);
                    % Distribution 2
                    distInfo2 = makedist(distributionLabel1,...
                        'Mu', Mu2, 'Sigma', Sigma2);
                    pdfCurve2 = pdf(distInfo2,obj.x);
                    % Distribution 3
                    distInfo3 = makedist(distributionLabel1,...
                        'Mu', Mu3, 'Sigma', Sigma3);
                    pdfCurve3 = pdf(distInfo3,obj.x);
                    % Mixture PDF Curve
                    obj.pdf_y = p(1)*pdfCurve1 +...
                        p(2)*pdfCurve2 + p(3)*pdfCurve3;
                    %\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
                    % generate random sample or actual pdf
                    if obj.randomVSactual == "random"
                        % mixture string array flag for misc_functions.mixSampling()
                        mixtureType = "three";
                        % generate n vector for mixture samplings
                        n = misc_functions.mixSampling(obj.Ns,p,mixtureType);
                        % generate random sample
                        rndData1 = random(distInfo1,1,n(1));
                        rndData2 = random(distInfo2,1,n(2));
                        rndData3 = random(distInfo3,1,n(3));
                        obj.rndData = [rndData1,rndData2,rndData3];
                    elseif obj.randomVSactual == "actual"
                        obj.actual_data = vertcat(obj.x,obj.pdf_y);
                    end
                    % CREATE DISTRIBUTION OBJECT --------------------------
                    % mixture string array flag for misc_functions.mixSampling()
                    mixtureType = "three";
                    % generate m vector for mixture samplings
                    m = misc_functions.mixSampling(obj.Ncdf,p,mixtureType);
                    % generate random sample to create distribution
                    % object
                    actData1 = random(distInfo1,1,m(1));
                    actData2 = random(distInfo2,1,m(2));
                    actData3 = random(distInfo3,1,m(3));
                    obj.actual_data = [actData1,actData2,actData3];
                    % generate numerical cdf: f=cdf, s=x-coordinates
                    [f,s] = ecdf(obj.actual_data);
                    f = f(1:2:end);
                    s = s(1:2:end);
                    % generate distribution object
                    obj.distInfo = ...
                        makedist('PiecewiseLinear','x',s,'Fx',f);
                    obj.cdf_y = cdf(obj.distInfo,obj.x);               
                    % create cdf/pdf from distribution object. 
                    % for debugging and vizualization.
                    if debug
                        if obj.randomVSactual == "actual"
                            xMix = linspace(obj.min_limit,...
                                obj.max_limit,obj.Npdf);
                            CDF = cdf(obj.distInfo,xMix);
                            % numerically differentiate
                            PDF = zeros(1,size(CDF(1:end-1),2));
                            for i = 1:size(CDF,2)-1
                                dx1 = (xMix(i+1)-xMix(i));
                                PDF(i) =  (CDF(i+1) - CDF(i))/dx1;
                            end
                            % smooth pdf obj.actual_data
                            smoo1 = smooth(xMix(1:end-1),PDF,obj.smooth);
                            % plot cdf,pdf,smoothed-pdf
                            figure('Name',['Debug: ',...
                                char(obj.dist_name)])
                            subplot(2,1,1)
                            hold on
                            plot(xMix,CDF,'-m')
                            plot(xMix(1:end-1),PDF,'-r')
                            plot(xMix(1:end-1),smoo1,'-b')
                            ylabel('$f(x) or F(x)$','Interpreter','latex')
                            xlabel('x','Interpreter','latex')
                            legend('cdf','pdf','smoothed-pdf')
                            
                            % plot histogram for random sample
                            subplot(2,1,2)
                            histogram(random(obj.distInfo,obj.Nhist,1),...
                                'Normalization','probability')
                            xlim([obj.min_limit,obj.max_limit])
                        end
                    end
                    %------------------------------------------------------
                case 'tLocationScale'
                    % t-Location Scale Case Statement
                    % Location parameter
                    Mu = 4;
                    % Scale parameter
                    Sigma = 0.05;
                    % Shape parameter
                    Nu = 1;
                    % Location parameter
                    % Delta = 3;
                    % PDF Curve \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
                    obj.distInfo = makedist(obj.dist_name,...
                        'Mu', Mu, 'Sigma', Sigma, 'Nu', Nu);
                    obj.pdf_y = pdf(obj.distInfo, obj.x);
                    obj.rndData = random(obj.distInfo,1,obj.Ns);
                    obj.actual_data = vertcat(obj.x,obj.pdf_y);
                    % CDF and Q Curve \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
                    obj.cdf_y = cdf(obj.distInfo,obj.x);
                    obj.true_quantile = icdf(obj.distInfo,obj.pVector);
                    %\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
                case 'Uniform'
                    % Uniform Case Statement
                    % Lower bound
                    Lower = 4;
                    % Upper Bound
                    Upper = 8;
                    % PDF Curve \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
                    obj.distInfo = makedist(obj.dist_name,...
                        'Lower', Lower, 'Upper', Upper);
                    obj.pdf_y = pdf(obj.distInfo,obj.x);
                    
                    obj.cdf_y = cdf(obj.distInfo,obj.x);
                    obj.true_quantile = icdf(obj.distInfo,obj.pVector);
                    obj.rndData = random(obj.distInfo,1,obj.Ns);
                    obj.actual_data = vertcat(obj.x,obj.pdf_y);
                    % CDF and Q Curve \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
                    obj.cdf_y = cdf(obj.distInfo,obj.x);
                    obj.true_quantile = icdf(obj.distInfo,obj.pVector);
                    %\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
                case 'Uniform-Mix' %------------------------------------------------- Mixture****
                    % Uniform Case Statement
                    % mixture weights
                    p1 = 0.1;
                    p2 = 0.6;
                    p3 = 1 - p1 - p2;
                    p = [p1,p2,p3];
                    % Lower bound
                    Lower1 = 1;
                    Lower2 = 3.5;
                    Lower3 = 7;
                    % Upper Bound
                    Upper1 = 2;
                    Upper2 = 5.5;
                    Upper3 = 9;
                    % PDF Curve \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
                    % Distribution 1
                    distributionLabel1 = 'Uniform';
                    distInfo1 = makedist(distributionLabel1,...
                        'Lower', Lower1, 'Upper', Upper1);
                    pdfCurve1 = pdf(distInfo1,obj.x);
                    % Distribution 2
                    distInfo2 = makedist(distributionLabel1,...
                        'Lower', Lower2, 'Upper', Upper2);
                    pdfCurve2 = pdf(distInfo2,obj.x);
                    % Distribution 3
                    distInfo3 = makedist(distributionLabel1,...
                        'Lower', Lower3, 'Upper', Upper3);
                    pdfCurve3 = pdf(distInfo3,obj.x);
                    % Mixture PDF Curve
                    obj.pdf_y = p(1)*pdfCurve1 +...
                        p(2)*pdfCurve2 + p(3)*pdfCurve3;
                    %\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
                    % generate random sample or actual pdf
                    if obj.randomVSactual == "random"
                        % mixture string array flag for misc_functions.mixSampling()
                        mixtureType = "three";
                        % generate n vector for mixture samplings
                        n = misc_functions.mixSampling(obj.Ns,p,mixtureType);
                        % generate random sample
                        rndData1 = random(distInfo1,1,n(1));
                        rndData2 = random(distInfo2,1,n(2));
                        rndData3 = random(distInfo3,1,n(3));
                        obj.rndData = [rndData1,rndData2,rndData3];
                    elseif obj.randomVSactual == "actual"
                        obj.actual_data = vertcat(obj.x,obj.pdf_y);
                    end
                    % CREATE DISTRIBUTION OBJECT --------------------------
                    % mixture string array flag for misc_functions.mixSampling()
                    mixtureType = "three";
                    % generate m vector for mixture samplings
                    m = misc_functions.mixSampling(obj.Ncdf,p,mixtureType);
                    % generate random sample to create distribution
                    % object
                    actData1 = random(distInfo1,1,m(1));
                    actData2 = random(distInfo2,1,m(2));
                    actData3 = random(distInfo3,1,m(3));
                    obj.actual_data = [actData1,actData2,actData3];
                    % generate numerical cdf: f=cdf, s=x-coordinates
                    [f,s] = ecdf(obj.actual_data);
                    f = f(1:2:end);
                    s = s(1:2:end);
                    % generate distribution object
                    obj.distInfo = ...
                        makedist('PiecewiseLinear','x',s,'Fx',f);
                    obj.cdf_y = cdf(obj.distInfo,obj.x);   
                    % create cdf/pdf from distribution object. 
                    % for debugging and vizualization.
                    if debug
                        if obj.randomVSactual == "actual"
                            xMix = linspace(obj.min_limit,...
                                obj.max_limit,obj.Npdf);
                            CDF = cdf(obj.distInfo,xMix);
                            % numerically differentiate
                            PDF = zeros(1,size(CDF(1:end-1),2));
                            for i = 1:size(CDF,2)-1
                                dx1 = (xMix(i+1)-xMix(i));
                                PDF(i) =  (CDF(i+1) - CDF(i))/dx1;
                            end
                            % smooth pdf obj.actual_data
                            smoo1 = smooth(xMix(1:end-1),PDF,obj.smooth);
                            % plot cdf,pdf,smoothed-pdf
                            figure('Name',['Debug: ',...
                                char(obj.dist_name)])
                            subplot(2,1,1)
                            hold on
                            plot(xMix,CDF,'-m')
                            plot(xMix(1:end-1),PDF,'-r')
                            plot(xMix(1:end-1),smoo1,'-b')
                            ylabel('$f(x) or F(x)$','Interpreter','latex')
                            xlabel('x','Interpreter','latex')
                            legend('cdf','pdf','smoothed-pdf')
                            % plot histogram for random sample
                            subplot(2,1,2)
                            histogram(random(obj.distInfo,obj.Nhist,1),...
                                'Normalization','probability')
                            xlim([obj.min_limit,obj.max_limit])
                        end
                    end
                    %------------------------------------------------------
                case 'Weibull'
                    % Weibull Case Statement
                    % Scale parameter
                    a = 1;
                    % Shape parameter
                    b = 2;
                    % PDF Curve \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
                    obj.distInfo = ...
                        makedist(obj.dist_name,'a', a, 'b', b);
                    obj.pdf_y = pdf(obj.distInfo,obj.x);
                    obj.rndData = random(obj.distInfo,1,obj.Ns);
                    obj.actual_data = vertcat(obj.x,obj.pdf_y);
                    % CDF and Q Curve \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
                    obj.cdf_y = cdf(obj.distInfo,obj.x);
                    obj.true_quantile = icdf(obj.distInfo,obj.pVector);
                    %\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
                otherwise
                    % Warning Statement
                    warning('No distribution was picked')
            end
            
            if obj.generate_data
                % Create obj.actual_data file
                if obj.randomVSactual == "random"
                    dataCreation(obj.rndData,obj.filename,obj.precision,1)
                elseif obj.randomVSactual == "actual"
                    dataCreation(obj.actual_data,obj.filename,obj.precision,1)
                end
                % Create folder for distribution obj.actual_data \\\\\\\\\\\\\\
                % Define folder name
                folderName = sprintf(['D_', char(obj.dist_name)]);
                % If folder already exist don't make it again
                if exist(folderName,'dir') == 0
                    mkdir(char(folderName))
                end
                %\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
                % Move datafile to folder
                if exist([char(obj.filename),'.txt'],'file') == 2
                    movefile([char(obj.filename),'.txt'] ,char(folderName));
                end
            end
        end
    end
end