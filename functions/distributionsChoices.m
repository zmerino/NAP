function pdfCurve = distributionsChoices(distributionName,x,fileName,randomVSactual,precision,Ns)
% Probability Distribution Data Generation function
% Created By: Zach D. Merino a MS candidate
% Updated: 3/12/19
%--------------------------------------------------------------------------
% This function generates the random sample for a distribution present in
% the list below as well as creates: sample data files, distribution data
% folders, and moves all files associated with a particular distribution to
% that folder.
%--------------------------------------------------------------------------
% There are case sensative operations used for distributionName
switch distributionName
    case 'Beta-a0p5-b1p5'
        %% Beta1 Case Statement
        
        % First shape parameter
        a = 0.5;
        % Second shape parameter
        b = 1.5;
        
        % PDF Curve \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
        probabilityDistribution = makedist('Beta','a',a,'b',b);
        pdfCurve = pdf(probabilityDistribution,x);
        %\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
        
        % generate random sample or actual pdf
        if randomVSactual == "random"
            rndData = random(probabilityDistribution,1,Ns);
        elseif randomVSactual == "actual"
            data = vertcat(x,pdfCurve);
        end
    case 'Beta-a2-b0p5'
        %% Beta2 Case Statement
        
        % First shape parameter
        a = 2;
        % Second shape parameter
        b = 0.5;
        
        % PDF Curve \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
        probabilityDistribution = makedist('Beta','a',a,'b',b);
        pdfCurve = pdf(probabilityDistribution,x);
        %\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
        
        % generate random sample or actual pdf
        if randomVSactual == "random"
            rndData = random(probabilityDistribution,1,Ns);
        elseif randomVSactual == "actual"
            data = vertcat(x,pdfCurve);
        end
    case 'Beta-a0p5-b0p5'
        %% Beta3 Case Statement
        
        % First shape parameter
        a = 0.5;
        % Second shape parameter
        b = 0.5;
        
        % PDF Curve \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
        probabilityDistribution = makedist('Beta','a',a,'b',b);
        pdfCurve = pdf(probabilityDistribution,x);
        %\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
        
        % generate random sample or actual pdf
        if randomVSactual == "random"
            rndData = random(probabilityDistribution,1,Ns);
        elseif randomVSactual == "actual"
            data = vertcat(x,pdfCurve);
        end
    case 'Bimodal-Normal'
        %% Normal Case Statement
        
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
        
        % PDF Curve \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
        % Distribution 1
        distributionLabel1 = 'Normal';
        probabilityDistribution1 = makedist(distributionLabel1,'Mu', Mu1, 'Sigma', Sigma1);
        pdfCurve1 = pdf(probabilityDistribution1,x);
        
        % Distribution 2
        probabilityDistribution2 = makedist(distributionLabel1,'Mu', Mu2, 'Sigma', Sigma2);
        pdfCurve2 = pdf(probabilityDistribution2,x);
        
        % Mixture PDF Curve
        pdfCurve = p(1)*pdfCurve1 + p(2)*pdfCurve2;
        data = vertcat(x,pdfCurve);
        %\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
        
        % generate random sample or actual pdf
        if randomVSactual == "random"
            % mixture string array flag for mixSampling()
            mixtureType = "two";
            % generate n vector for mixture samplings
            n = mixSampling(Ns,p,mixtureType);
            % generate random sample
            rndData1 = random(probabilityDistribution1,1,n(1));
            rndData2 = random(probabilityDistribution2,1,n(2));
            rndData = [rndData1,rndData2];
        elseif randomVSactual == "actual"
            data = vertcat(x,pdfCurve);
        end
    case 'Binomial'
        %% Binomial Case Statement
        
        % Number of trials
        n = 2000;
        % Porbability of success for each trial
        p = 0.2;
        
        % PDF Curve \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
        probabilityDistribution = makedist(distributionName,'n',n,'p',p);
        pdfCurve = binopdf(x,n,p);
        %\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
        
        % generate random sample or actual pdf
        if randomVSactual == "random"
            rndData = binornd(Ns,p);
        elseif randomVSactual == "actual"
            data = vertcat(x,pdfCurve);
        end
    case 'BirnbaumSaunders'
        %% BirnbaumSaunders Case Statement
        
        % Scale parameter
        Beta = 1.5;
        % Shape parameter
        Gamma = 0.5;
        
        % PDF Curve \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
        probabilityDistribution = makedist(distributionName,'Beta',Beta,'Gamma',Gamma);
        pdfCurve = pdf(probabilityDistribution,x);
        %\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
        
        % generate random sample or actual pdf
        if randomVSactual == "random"
            rndData = random(probabilityDistribution,1,Ns);
        elseif randomVSactual == "actual"
            data = vertcat(x,pdfCurve);
        end
        
    case 'BirnbaumSaunders-Stable'
        %% BirnbaumSaunders Case Statement
        
        % mixture weights
        p1 = 0.35;
        p2 = 1 - p1;
        p = [p1,p2];
        
        % BirnbaumSaunders distribution -----------------------------------
        % Scale parameter
        Beta = 1.5;
        % Shape parameter
        Gamma = 0.5;
        
        % Stable distribution ---------------------------------------------
        % First shape parameter
        Alpha1 = 0.5;
        % Second shape parameter: -1 <= Beta <= 1
        Beta1 = 0.05;
        % Scale parameter
        Gam1 = 1;
        % Location parameter
        Delta1 = 7;
        
        % PDF Curve \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
        % BirnbaumSaunders distribution
        distributionLabel1 = 'BirnbaumSaunders';
        probabilityDistribution1 = makedist(distributionLabel1,'Beta',Beta,'Gamma',Gamma);
        pdfCurve1 = pdf(probabilityDistribution1,x);
        
        % Stable distribution
        distributionLabel2 = 'Stable';
        probabilityDistribution2 = makedist(distributionLabel2,'Alpha', Alpha1,...
            'Beta', Beta1, 'Gam', Gam1, 'Delta', Delta1);
        pdfCurve2 = pdf(probabilityDistribution2, x);
        
        % Mixture PDF Curve
        pdfCurve = p(1)*pdfCurve1 + p(2)*pdfCurve2;
        data = vertcat(x,pdfCurve);
        %\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
        
        % generate random sample or actual pdf
        if randomVSactual == "random"
            % mixture string array flag for mixSampling()
            mixtureType = "two";
            % generate n vector for mixture samplings
            n = mixSampling(Ns,p,mixtureType);
            % generate random sample
            rndData1 = random(probabilityDistribution1,1,n(1));
            rndData2 = random(probabilityDistribution2,1,n(2));
            rndData = [rndData1,rndData2];
        elseif randomVSactual == "actual"
            data = vertcat(x,pdfCurve);
        end
        
    case 'Burr'
        %% Burr Case Statement
        
        % Scale parameter
        Alpha = 1;
        % Shape parameter one
        c = 2;
        % Shape parameter one two
        k = 2;
        
        % PDF Curve \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
        probabilityDistribution = makedist(distributionName,'Alpha',Alpha,'c',c,'k',k);
        pdfCurve = pdf(probabilityDistribution,x);
        %\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
        
        % generate random sample or actual pdf
        if randomVSactual == "random"
            rndData = random(probabilityDistribution,1,Ns);
        elseif randomVSactual == "actual"
            data = vertcat(x,pdfCurve);
        end
    case 'Chisquare'
        %% Chisquare Case Statement
        
        % Degrees of freedom
        Nu = 4;
        
        % PDF Curve \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
        pdfCurve = chi2pdf(x,Nu);
        %\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
        
        % generate random sample or actual pdf
        if randomVSactual == "random"
            rndData = chi2rnd(Nu,1,Ns);
        elseif randomVSactual == "actual"
            data = vertcat(x,pdfCurve);
        end
    case 'Exponential'
        %% Exponential Case Statement
        
        % Mean
        Mu = 1;
        
        % PDF Curve \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
        probabilityDistribution = makedist(distributionName,'Mu',Mu);
        pdfCurve = pdf(probabilityDistribution,x);
        %\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
        
        % generate random sample or actual pdf
        if randomVSactual == "random"
            rndData = random(probabilityDistribution,1,Ns);
        elseif randomVSactual == "actual"
            data = vertcat(x,pdfCurve);
        end
    case 'Extreme-Value'
        %% Extreme Value Case Statement
        
        % Location parameter
        Mu = 1;
        % Scale parameter
        Sigma = 2;
        
        % PDF Curve \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
        distributionLabel = 'Extreme Value';
        probabilityDistribution = makedist(distributionLabel,'Mu',Mu, 'Sigma', Sigma);
        pdfCurve = pdf(probabilityDistribution,x);
        %\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
        
        % generate random sample or actual pdf
        if randomVSactual == "random"
            rndData = random(probabilityDistribution,1,Ns);
        elseif randomVSactual == "actual"
            data = vertcat(x,pdfCurve);
        end
    case 'Gamma'
        %% Gamma Case Statement
        
        % Shape parameter
        a = 2;
        % Scale parameter
        b = 2;
        
        % PDF Curve \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
        probabilityDistribution = makedist(distributionName,'a',a,'b',b);
        pdfCurve = pdf(probabilityDistribution,x);
        %\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
        
        % generate random sample or actual pdf
        if randomVSactual == "random"
            rndData = random(probabilityDistribution,1,Ns);
        elseif randomVSactual == "actual"
            data = vertcat(x,pdfCurve);
        end
    case 'Generalized-Extreme-Value'
        %% Generalized Extreme Value Value Statement
        
        % Shape parameter
        k = 1;
        % Scale parameter
        Sigma = 2;
        % Location parameter
        Mu = 2;
        
        % PDF Curve \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
        distributionLabel = 'Generalized Extreme Value';
        probabilityDistribution = makedist(distributionLabel,'k',k, 'Sigma', Sigma,'Mu',Mu);
        pdfCurve = pdf(probabilityDistribution,x);
        %\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
        
        % generate random sample or actual pdf
        if randomVSactual == "random"
            rndData = random(probabilityDistribution,1,Ns);
        elseif randomVSactual == "actual"
            data = vertcat(x,pdfCurve);
        end
    case 'Generalized-Pareto'
        %% Generalized Pareto Value Value Case Statement
        
        % Tail MemTracker (shape) parameter
        k = 2;
        % Scale parameter
        Sigma = 1;
        % Threshold (location) parameter
        theta = 0;
        
        % PDF Curve \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
        distributionLabel = 'Generalized Pareto';
        probabilityDistribution = makedist(distributionLabel,'k',k, 'Sigma', Sigma,'Theta', theta);
        pdfCurve = pdf(probabilityDistribution,x);
        % Mu parameter is not recognized
        %\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
        
        % generate random sample or actual pdf
        if randomVSactual == "random"
            rndData = random(probabilityDistribution,1,Ns);
        elseif randomVSactual == "actual"
            data = vertcat(x,pdfCurve);
        end
    case 'HalfNormal'
        %% Half Normal Value Case Statement
        
        % Location parameter
        Mu = 0;
        % Scale parameter
        Sigma = 1;
        
        % PDF Curve \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
        probabilityDistribution = makedist(distributionName,'Mu', Mu, 'Sigma', Sigma);
        pdfCurve = pdf(distributionName,x);
        %\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
        
        % generate random sample or actual pdf
        if randomVSactual == "random"
            rndData = random(probabilityDistribution,1,Ns);
        elseif randomVSactual == "actual"
            data = vertcat(x,pdfCurve);
        end
    case 'InverseGaussian'
        %% Inverse Gaussian Case Statement
        
        % Scale parameter
        Mu = 1;
        % Shape parameter
        Lambda = 1;
        
        % PDF Curve \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
        probabilityDistribution = makedist(distributionName,'mu', Mu, 'lambda', Lambda);
        pdfCurve = pdf(probabilityDistribution,x);
        %\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
        
        % generate random sample or actual pdf
        if randomVSactual == "random"
            rndData = random(probabilityDistribution,1,Ns);
        elseif randomVSactual == "actual"
            data = vertcat(x,pdfCurve);
        end
    case 'Mix-Beta-Stable-1'
        % mixture model for stable and beta distributions
        
        % mixture weights
        p1 = 0.25;
        p2 = 1 - p1;
        p = [p1,p2];
        
        % Stable distribution ---------------------------------------------
        % First shape parameter
        Alpha1 = 0.3;
        %{
        Alpha = 0.4;
        Alpha = 0.35
        Alpha = 0.5  %<---- use to start
        Alpha = 0.2; %<---- very hard to estimate
        %}
        % Second shape parameter: -1 <= Beta <= 1
        Beta1 = 0.05;
        %{
        Beta = 0.9;
        Beta = 1;
        Beta = .05;
        Beta = .05;
        %}
        Gam1 = 1;
        % Location parameter
        Delta1 = 0.5;
        
        % Beta distribution -----------------------------------------------
        % First shape parameter
        a = 0.5;
        % Second shape parameter
        b = 0.5;
        
        % PDF Curve \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
        % Distribution 1
        distributionLabel1 = 'Beta';
        probabilityDistribution1 = makedist(distributionLabel1,'a',a,'b',b);
        pdfCurve1 = pdf(probabilityDistribution1, x);
        
        % Distribution 2
        distributionLabel2 = 'Stable';
        probabilityDistribution2 = makedist(distributionLabel2,'Alpha', Alpha1,...
            'Beta', Beta1, 'Gam', Gam1, 'Delta', Delta1);
        pdfCurve2 = pdf(probabilityDistribution2, x);
        
        % Mixture PDF Curve
        pdfCurve = p(1)*pdfCurve1 + p(2)*pdfCurve2;
        %\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
        
        % generate random sample or actual pdf
        if randomVSactual == "random"
            % mixture string array flag for mixSampling()
            mixtureType = "two";
            % generate n vector for mixture samplings
            n = mixSampling(Ns,p,mixtureType);
            % generate random sample
            rndData1 = random(probabilityDistribution1,1,n(1));
            rndData2 = random(probabilityDistribution2,1,n(2));
            rndData = [rndData1,rndData2];
        elseif randomVSactual == "actual"
            data = vertcat(x,pdfCurve);
        end
        
    case 'Normal'
        %% Normal Case Statement
        
        % Mean
        Mu = 5;
        % Standard deviation
        Sigma = 1;
        
        % PDF Curve \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
        probabilityDistribution = makedist(distributionName,'Mu', Mu, 'Sigma', Sigma);
        pdfCurve = pdf(probabilityDistribution,x);
        %\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
        
        % generate random sample or actual pdf
        if randomVSactual == "random"
            rndData = random(probabilityDistribution,1,Ns);
        elseif randomVSactual == "actual"
            data = vertcat(x,pdfCurve);
        end
        
    case 'Normal-Contaminated'
        %% Normal Case Statement
        
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
        
        % PDF Curve \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
        % Distribution 1
        distributionLabel1 = 'Normal';
        distInfo1 = makedist(distributionLabel1,'Mu', Mu1, 'Sigma', Sigma1);
        pdfCurve1 = pdf(distInfo1,x);
        
        % Distribution 2
        distInfo2 = makedist(distributionLabel1,'Mu', Mu2, 'Sigma', Sigma2);
        pdfCurve2 = pdf(distInfo2,x);
        
        % Mixture PDF Curve
        pdfCurve = p(1)*pdfCurve1 + p(2)*pdfCurve2;
        data = vertcat(x,pdfCurve);
        %\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
        
        % generate random sample or actual pdf
        if randomVSactual == "random"
            % mixture string array flag for mixSampling()
            mixtureType = "two";
            % generate n vector for mixture samplings
            n = mixSampling(Ns,p,mixtureType);
            % generate random sample
            rndData1 = random(distInfo1,1,n(1));
            rndData2 = random(distInfo2,1,n(2));
            rndData = [rndData1,rndData2];
        elseif randomVSactual == "actual"
            data = vertcat(x,pdfCurve);
        end
    case 'Poisson'
        %% Poisson Case Statement
        
        % Mean
        Lambda = 0.001;
        
        % PDF Curve \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
        probabilityDistribution = makedist(distributionName,'Lambda', Lambda);
        % Graphing issue or Lambda parameter incorrectly picked?
        pdfCurve = poisspdf(x,Lambda);
        %\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
        
        % generate random sample or actual pdf
        if randomVSactual == "random"
            rndData = random(probabilityDistribution,1,Ns);
        elseif randomVSactual == "actual"
            data = vertcat(x,pdfCurve);
        end
        
    case 'Square-periodic'
        %% Uniform Case Statement
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
        
        
        % PDF Curve \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
        distributionLabel1 = 'Uniform';
        % Distribution 1
        distInfo1 = makedist(distributionLabel1,'Lower', Lower1, 'Upper', Upper1);
        pdfCurve1 = pdf(distInfo1,x);
        
        % Distribution 2
        distInfo2 = makedist(distributionLabel1,'Lower', Lower2, 'Upper', Upper2);
        pdfCurve2 = pdf(distInfo2,x);
        
        % Distribution 3
        distInfo3 = makedist(distributionLabel1,'Lower', Lower3, 'Upper', Upper3);
        pdfCurve3 = pdf(distInfo3,x);
        
        % Distribution 4
        distInfo4 = makedist(distributionLabel1,'Lower', Lower4, 'Upper', Upper4);
        pdfCurve4 = pdf(distInfo4,x);
        
        % Distribution 5
        distInfo5 = makedist(distributionLabel1,'Lower', Lower5, 'Upper', Upper5);
        pdfCurve5 = pdf(distInfo5,x);
        
        % Distribution 6
        distInfo6 = makedist(distributionLabel1,'Lower', Lower6, 'Upper', Upper6);
        pdfCurve6 = pdf(distInfo6,x);
        
        % Mixture PDF Curve
        pdfCurve = p(1)*pdfCurve1 + p(2)*pdfCurve2 + p(3)*pdfCurve3 + p(4)*pdfCurve4 + p(5)*pdfCurve5 + p(6)*pdfCurve6;
        %\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
        
        % generate random sample or actual pdf
        if randomVSactual == "random"
            % mixture string array flag for mixSampling()
            mixtureType = "six";
            % generate n vector for mixture samplings
            n = mixSampling(Ns,p,mixtureType);
            % generate random sample
            rndData1 = random(distInfo1,1,n(1));
            rndData2 = random(distInfo2,1,n(2));
            rndData3 = random(distInfo3,1,n(3));
            rndData4 = random(distInfo4,1,n(4));
            rndData5 = random(distInfo5,1,n(5));
            rndData6 = random(distInfo6,1,n(6));
            rndData = [rndData1,rndData2,rndData3,rndData4,rndData5,rndData6];
        elseif randomVSactual == "actual"
            data = vertcat(x,pdfCurve);
        end
        
    case 'Stable'
        %% Stable Case Statement
        
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
        
        % PDF Curve \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
        probabilityDistribution = makedist(distributionName,'Alpha', Alpha,...
            'Beta', Beta, 'Gam', Gam, 'Delta', Delta);
        pdfCurve = pdf(probabilityDistribution, x);
        %\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
        
        % generate random sample or actual pdf
        if randomVSactual == "random"
            rndData = random(probabilityDistribution,1,Ns);
        elseif randomVSactual == "actual"
            data = vertcat(x,pdfCurve);
        end
    case 'Stable1'
        %% Stable Case Statement
        
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
        
        % PDF Curve \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
        distributionLabel = 'Stable';
        
        distInfo = makedist(distributionLabel,'Alpha', Alpha,...
            'Beta', Beta, 'Gam', Gam, 'Delta', Delta);
        pdfCurve = pdf(distInfo, x);
        %\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
        
        % generate random sample or actual pdf
        if randomVSactual == "random"
            rndData = random(distInfo,1,Ns);
        elseif randomVSactual == "actual"
            data = vertcat(x,pdfCurve);
        end
        
    case 'Stable2'
        %% mixture model for 2 stable distributions
        
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
        
        % PDF Curve \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
        distributionLabel = 'Stable';
        % stable 1
        distInfo1 = makedist(distributionLabel,'Alpha', Alpha1,...
            'Beta', Beta1, 'Gam', Gam1, 'Delta', Delta1);
        pdfCurve1 = pdf(distInfo1, x);
        
        % stable 2
        distInfo2 = makedist(distributionLabel,'Alpha', Alpha2,...
            'Beta', Beta2, 'Gam', Gam2, 'Delta', Delta2);
        pdfCurve2 = pdf(distInfo2, x);
        
        % Mixture PDF Curve
        pdfCurve = p(1)*pdfCurve1 + p(2)*pdfCurve2;
        %\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
        
        % generate random sample or actual pdf
        if randomVSactual == "random"
            % mixture string array flag for mixSampling()
            mixtureType = "two";
            % generate n vector for mixture samplings
            n = mixSampling(Ns,p,mixtureType);
            % generate random sample
            rndData1 = random(distInfo1,1,n(1));
            rndData2 = random(distInfo2,1,n(2));
            rndData = [rndData1,rndData2];
        elseif randomVSactual == "actual"
            data = vertcat(x,pdfCurve);
        end
        
    case 'Stable2'
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
        
        % PDF Curve \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
        distributionLabel = 'Stable';
        % stable 1
        probabilityDistribution1 = makedist(distributionLabel,'Alpha', Alpha1,...
            'Beta', Beta1, 'Gam', Gam1, 'Delta', Delta1);
        pdfCurve1 = pdf(probabilityDistribution1, x);
        
        % stable 2
        probabilityDistribution2 = makedist(distributionLabel,'Alpha', Alpha2,...
            'Beta', Beta2, 'Gam', Gam2, 'Delta', Delta2);
        pdfCurve2 = pdf(probabilityDistribution2, x);
        
        % Mixture PDF Curve
        pdfCurve = p(1)*pdfCurve1 + p(2)*pdfCurve2;
        %\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
        
        % generate random sample or actual pdf
        if randomVSactual == "random"
            % mixture string array flag for mixSampling()
            mixtureType = "two";
            % generate n vector for mixture samplings
            n = mixSampling(Ns,p,mixtureType);
            % generate random sample
            rndData1 = random(probabilityDistribution1,1,n(1));
            rndData2 = random(probabilityDistribution2,1,n(2));
            rndData = [rndData1,rndData2];
        elseif randomVSactual == "actual"
            data = vertcat(x,pdfCurve);
        end
    case 'Stable3'
        % mixture model for 2 stable distributions
        
        % mixture weights
        p1 = 0.25;
        p2 = 0.5;
        p3 = 1 - p1 - p2;
        p = [p1,p2,p3];
        
        % Stable distributions --------------------------------------------
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
        
        % PDF Curve \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
        distributionLabel = 'Stable';
        % stable 1
        probabilityDistribution1 = makedist(distributionLabel,'Alpha', Alpha1,...
            'Beta', Beta1, 'Gam', Gam1, 'Delta', Delta1);
        pdfCurve1 = pdf(probabilityDistribution1, x);
        
        % stable 2
        probabilityDistribution2 = makedist(distributionLabel,'Alpha', Alpha2,...
            'Beta', Beta2, 'Gam', Gam2, 'Delta', Delta2);
        pdfCurve2 = pdf(probabilityDistribution2, x);
        
        % stable 3
        probabilityDistribution3 = makedist(distributionLabel,'Alpha', Alpha3,...
            'Beta', Beta3, 'Gam', Gam3, 'Delta', Delta3);
        pdfCurve3 = pdf(probabilityDistribution3, x);
        
        % Mixture PDF Curve
        pdfCurve = p(1)*pdfCurve1 + p(2)*pdfCurve2 + p(3)*pdfCurve3;
        %\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
        
        % generate random sample or actual pdf
        if randomVSactual == "random"
            % mixture string array flag for mixSampling()
            mixtureType = "three";
            % generate n vector for mixture samplings
            n = mixSampling(Ns,p,mixtureType);
            % generate random sample
            rndData1 = random(probabilityDistribution1,1,n(1));
            rndData2 = random(probabilityDistribution2,1,n(2));
            rndData3 = random(probabilityDistribution3,1,n(3));
            rndData = [rndData1,rndData2,rndData3];
        elseif randomVSactual == "actual"
            data = vertcat(x,pdfCurve);
        end
    case 'Trimodal-Normal'
        %% Normal Case Statement
        
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
        
        % PDF Curve \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
        % Distribution 1
        distributionLabel1 = 'Normal';
        probabilityDistribution1 = makedist(distributionLabel1,'Mu', Mu1, 'Sigma', Sigma1);
        pdfCurve1 = pdf(probabilityDistribution1,x);
        
        % Distribution 2
        probabilityDistribution2 = makedist(distributionLabel1,'Mu', Mu2, 'Sigma', Sigma2);
        pdfCurve2 = pdf(probabilityDistribution2,x);
        
        % Distribution 3
        probabilityDistribution3 = makedist(distributionLabel1,'Mu', Mu3, 'Sigma', Sigma3);
        pdfCurve3 = pdf(probabilityDistribution3,x);
        
        % Mixture PDF Curve
        pdfCurve = p(1)*pdfCurve1 + p(2)*pdfCurve2 + p(3)*pdfCurve3;
        %\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
        
        % generate random sample or actual pdf
        if randomVSactual == "random"
            % mixture string array flag for mixSampling()
            mixtureType = "three";
            % generate n vector for mixture samplings
            n = mixSampling(Ns,p,mixtureType);
            % generate random sample
            rndData1 = random(probabilityDistribution1,1,n(1));
            rndData2 = random(probabilityDistribution2,1,n(2));
            rndData3 = random(probabilityDistribution3,1,n(3));
            rndData = [rndData1,rndData2,rndData3];
        elseif randomVSactual == "actual"
            data = vertcat(x,pdfCurve);
        end
    case 'tLocationScale'
        %% t-Location Scale Case Statement
        
        % Location parameter
        Mu = 0.5;
        % Scale parameter
        Sigma = .05;
        % Shape parameter
        Nu = 1;
        % Location parameter
        % Delta = 3;
        
        % PDF Curve \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
        probabilityDistribution = makedist(distributionName,'Mu', Mu, 'Sigma', Sigma, 'Nu', Nu);
        pdfCurve = pdf(probabilityDistribution, x);
        %\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
        
        % generate random sample or actual pdf
        if randomVSactual == "random"
            rndData = random(probabilityDistribution,1,Ns);
        elseif randomVSactual == "actual"
            data = vertcat(x,pdfCurve);
        end
    case 'Uniform'
        %% Uniform Case Statement
        
        % Lower bound
        Lower = 4;
        % Upper Bound
        Upper = 8;
        
        % PDF Curve \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
        probabilityDistribution = makedist(distributionName,'Lower', Lower, 'Upper', Upper);
        pdfCurve = pdf(probabilityDistribution, x);
        %\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
        
        % generate random sample or actual pdf
        if randomVSactual == "random"
            rndData = random(probabilityDistribution,1,Ns);
        elseif randomVSactual == "actual"
            data = vertcat(x,pdfCurve);
        end
    case 'Uniform-Mix'
        %% Uniform Case Statement
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
        
        
        % PDF Curve \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
        % Distribution 1
        distributionLabel1 = 'Uniform';
        distInfo1 = makedist(distributionLabel1,'Lower', Lower1, 'Upper', Upper1);
        pdfCurve1 = pdf(distInfo1,x);
        
        % Distribution 2
        distInfo2 = makedist(distributionLabel1,'Lower', Lower2, 'Upper', Upper2);
        pdfCurve2 = pdf(distInfo2,x);
        
        % Distribution 3
        distInfo3 = makedist(distributionLabel1,'Lower', Lower3, 'Upper', Upper3);
        pdfCurve3 = pdf(distInfo3,x);
        
        % Mixture PDF Curve
        pdfCurve = p(1)*pdfCurve1 + p(2)*pdfCurve2 + p(3)*pdfCurve3;
        %\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
        
        % generate random sample or actual pdf
        if randomVSactual == "random"
            % mixture string array flag for mixSampling()
            mixtureType = "three";
            % generate n vector for mixture samplings
            n = mixSampling(Ns,p,mixtureType);
            % generate random sample
            rndData1 = random(distInfo1,1,n(1));
            rndData2 = random(distInfo2,1,n(2));
            rndData3 = random(distInfo3,1,n(3));
            rndData = [rndData1,rndData2,rndData3];
        elseif randomVSactual == "actual"
            data = vertcat(x,pdfCurve);
        end
        
    case 'Weibull'
        %% Weibull Case Statement
        
        % Scale parameter
        a = 1;
        % Shape parameter
        b = 2;
        
        % PDF Curve \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
        probabilityDistribution = makedist(distributionName,'a', a, 'b', b);
        pdfCurve = pdf(probabilityDistribution, x);
        %\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
        
        % generate random sample or actual pdf
        if randomVSactual == "random"
            rndData = random(probabilityDistribution,1,Ns);
        elseif randomVSactual == "actual"
            data = vertcat(x,pdfCurve);
        end
    otherwise
        %% Warning Statement
        warning('No distribution was picked')
end

% Create data file
if randomVSactual == "random"
    dataCreation(rndData,fileName,precision,1)
elseif randomVSactual == "actual"
    dataCreation(data,fileName,precision,1)
end

% Create folder for distribution data \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
% Define folder name
folderName = sprintf(['D_', char(distributionName)]);
% If folder already exist don't make it again
if exist(folderName,'dir') == 0
    mkdir(char(folderName))
end
%\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

% Move datafile to folder


if exist([char(fileName),'.txt'],'file') == 2 && ~contains(char(fileName),'stitchPDF')
    movefile([char(fileName),'.txt'] ,char(folderName));
end

end