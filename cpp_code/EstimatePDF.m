% EstimatePDF computes a probability density estimate for a one-dimensional
% data sample based on a nonparametric maximum entropy method.  For additional 
% usage information see readme.txt file with this installation.
% For conceptual details of the method and algorithm, see:
%
% Farmer, Jenny and Donald Jacobs (2018). 
% "High throughput nonparametric probability density estimation." PLoS One 13(5): e0196937.
%
%
% [FAILED, XI, F, CDF, SQR, LAGRANGE, SCORE, CONFIDENCE, SURD] = EstimatePDF(X) 
% Computes the density estimate of data in sample X with default settings.
% Outputs returned are as follows:
%
%   FAILED              non-zero if solution not found
%   XI                  relative spacing along x-axis for density estimate
%   F                   probability density function for estimate (PDF)
%   CDF                 cummulative density function for estimate
%   SQR                 scaled quantile residual for sample data
%   LAGRANGE            Lagrange multipliers estimated to construct PDF
%   SCORE               Value of final score for estimate
%   CONFIDENCE          Confidence threshold for estimate
%   SURD                sampled uniform random data
%
% [...] = EstimatePDF(X, PARM) Computes the density estimate of data in 
% sample X with parameters in a structure defined in PARM.  Parameter options 
% are listed below.  
% 
% Parameter name        Default value               Description
%
% lowBound              calculated                  set fixed lower bound
% highBound             calculated                  set fixed upper bound
% integrationPoints 	max(1500, 200/n + 200)      data resolution; n=sample size
% LagrangeMin           1                           minimum number of weighted expansions  
% LagrangeMax           200                         maximum number of weighted expansions 		
% SURDtarget            40                          target confidence threshold 
% SURDmin               5                           minimum conficence threshold accepted
% SURDmax               100                         maximum conficence threshold accepted
% scoreType             ‘QZ'                        Change to ‘LL’ for log likelihood scoring
% debug                 false                       detailed output to console
% partition             1025                        Initial data partition for scoring, set to zero for no partitioning
% outlierCutoff         7                           Set to zero to disable outlier detection
% adaptiveDx            true                        Set to false for uniformly spaced numerical solution
% 		
%
% Example:  Plot estimated density for a Normal distribution using default
% settings:
%   [failed, y, pdf] = EstimatePDF(randn(1000, 1));
%   if ~failed
%       plot(y, pdf);
%   end
%
% Example:  Plot estimated cummulative density for a uniform distribution defined on
% the interval (0 1):
%   parms.lowBound = 0;
%   parms.highBound = 1;
%   [failed, y, pdf, cdf] = EstimatePDF(rand(1000, 1), parms);
%   if ~failed
%       plot(y, cdf);
%   end
%
% Example:  Plot estimated  density for an exponential distribution defined on
% the interval (0 inf), requiring at least 2 Langrange multipliers:
%   parms.lowBound = 0;
%   parms.LagrangeMin = 2;
%   [failed, y, pdf] = EstimatePDF(exprnd(1, 1000, 1), parms);
%   if ~failed
%       plot(y, pdf);
%   end





