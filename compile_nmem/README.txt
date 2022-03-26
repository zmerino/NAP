================================================================
PDFAnalyze, Version 1.0, April 2020
Jenny Farmer jfarmer@carolina.rr.com
Donald Jacobs djacobs1@uncc.edu
University of North Carolina at Charlotte



================================================================
GENERAL INFORMATION
================================================================

PDFAnalyze computes a probability density estimate for a one-dimensional data sample and produces optional plots for analysis.  
It is based on a nonparametric probability density estimation method called PDFEstimator.

Please cite this publication if you use this code for your research:

https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0196937



=================================================================
INSTALLATION FOR MATLAB (2018 or greater)
=================================================================

Installation Steps


1. Prior to installing the PDFAnalyze, the MingGW C/C++ compiler for Windows must be installed as a MATLAB Add-on.  
To install, select [Add-Ons/Get Add-Ons] from the HOME menu within MATLAB and search for ‘MinGW’.  Select and install MinGW-w64.

2. Copy all source files into a single folder

3. Run the CompilePDF.m script in MATLAB to create a MATLAB Executable (mex) 

4. (optional) Verify installation by running example.m script in MATLAB 



The PDFAnalyze package consists of the following files:

PDFAnalyze.m
PlotBeta.m
EstimatePDF.m
FigureSettings.m
GetTargets.m
example.m; 
CompilePDF.m; 

EstimatePDF.cpp; EstimatePDF.h 
callPDF.cpp;  callPDF.h
ChebyShev.cpp; ChebyShev.h
InputData.cpp; InputData.h
InputParameters.cpp; InputParameters.h
MinimizeScore.cpp; MinimizeScore.h
Score.cpp; Score.h
ScoreQZ.cpp; ScoreQZ.h
ScoreLL.cpp; ScoreLL.h
WriteResults.cpp; WriteResults.h
OutputControl.cpp; OutputControl.h

README.txt



=================================================================
USAGE
=================================================================


[F, XI] = PDFAnalyze(X) Computes the density estimate of data in sample X.  F contains the density estimate at points XI.  
	The number of points and the relative spacing is determined automatically from the features of the data sample


[F, XI, CDF, SQR] = PDFAnalyze(...) also returns the cumulative density and the scaled quantile residual for each sample data point.


PDFAnalyze(...) with no output arguments produces a plot of the density estimate.


[...] = PDFAnalyze(..., 'param1', 'val1, 'param2', 'val2', ...) specifies parameter name/value pairs to control the density estimation.  

Valid parameters are as follows:

	Parameter			Value
	'PlotType'			Produces any combination of three plot types:
				'pdf'		probability density function for the
                	           		'sqr'		scaled quantile residual
				'combined'		pdf and sqr plotted on one figure

				Multiple plot types occur with multiple name/value pairs specified

	'EstimationType'		The default estimation method is PDFEstimate.

				Additional KDE methods are available:
				'kde1'		built-in MATLAB function ksdensity
				'kde2'		Zdravko Botev (2020). Kernel Density Estimator 
						(https://www.mathworks.com/matlabcentral/fileexchange/14034-kernel-density-est
						MATLAB Central File Exchange. Retrieved March 17, 2020.

	'distribution'		A two column matrix, [F, XI], representing a distribution to plot on the same figure as the estimate for use with 'pdf' plot type.  
				Useful for comparison to a known distribution.



=================================================================
EXAMPLES
=================================================================

Example 1: Plot the estimate of random sample for a Normal distribution along with the true Normal distribution:

	data = randn(1000, 1);
	x = min(data):0.1:max(data);
	f = normpdf(x);
	d = [x(:), f(:)];
	PDFAnalyze(data, 'distribution', d);


Example 2: Plot the scaled quantile residual (SQR) for an estimate of the Normal distribution, showing confidence thresholds and uncertainty estimates:

	PDFAnalyze(randn(10000, 1), 'PlotType', 'sqr');



=================================================================
PDFEstimator Usage
=================================================================

The PDFEstimator is invoked from within PDFAnalyze and can be customized through a collection of advanced input and output options.



Usage

[failed, y, pdf, cdf, sqr, lagrange, score, confidence, SURD] = EstimatePDF(data, parameters)

data		(required) a single vector of random sample data.
parameters	(optional) a MATLAB structure of options listed below



Optional Input Parameters

Name			Default Value
parameters.SURDtarget	[40]      
parameters.SURDmin		[5]
parameters.SURDmax	[100] 
parameters.LagrangeMin	[1]
parameters.LagrangeMax	[200]
parameters.lowBound		[calculated]
parameters.highBound	[calculated]      
parameters.integrationPoints	[calculated]
parameters.debug		[false]
parameters.partition		[1025]      
parameters.scoreType	['QZ']
parameters.outlierCutoff	[7]
parameters.adaptiveDx	[true]



Output Parameters

failed		non-zero if a solution was not found
y		range of values in PDF (independent variable)	
pdf		Probability Density Function (PDF)
cdf		Cummulative Denstiy Function
sqr		Scaled Quantile Residual
lagrange		Lagrange coefficients
score		Value returned by the score-type selected
confidence		SURD threshold achieved 
SURD		Sample Uniform Random Data            





NOTES

The following section includes a few brief notes concerning more advanced input and output options available, and how they may affect performance of the estimation.  
For a greater understanding of the methodology used, please see the publication referenced in the GENERAL INFORMATION section.


1. SURD Threshold Targets

Sample Uniform Random Data (SURD) loosely correlates with the strength of the solution, with higher thresholds indicating more probably solutions for the PDF.  


2. Scaled Quantile Residual

The equation for Scaled Quantile Residual (SQR) is given by SQR = sqrt(N+2)*(u - uniform-u) where N is the number of data samples. 
SQR plots are very useful as a diagnostic measure because they are sample size invariant and have universal characteristics independent of the true PDF.


3. Lagrange Coefficients

Each Lagrange multiplier returned as output is an expansion coefficient in the series of orthogonal functions within an exponential. 
The more complex the shape of the distribution, the more Lagrange multipliers are required to accurately define the PDF. 


4. Greater accuracy in numerical integration can be controlled

Increasing the number of integration points will improve the resolution of the PDF, but could increase runtime.  
Decreasing the integration points is not recommended, as it may produce poor solutions.  


5. Failed solutions 

Two safety measures are implemented to prevent the program from continuing an unreasonably long time without finding a solution. 

i) If progress stalls and the score is not improving significantly after many attempts, or 
ii) If the maximum number of Lagrange multipliers has been reached. 

If the maximum number of Lagrange multipliers is reached, this indicates that the solution is likely not yet converged. The user can increase the maximum. 
However, the default maximum of 200 is set to prevent cases that may never converge. 


6. Parametric maximum entropy method can be used with this program.

If a user desires an exact number of Lagrange multipliers, the  minLagrange and maxLagrange parameter options can be set to equal values. 
For example, if the user knows the distribution is a Gaussian, then the user could set both the minimum and maximum Lagrange mutlipliers to 3.  
In this case, the output will be equivalent to a parametric maximum entropy method, where the number of Lagrange multipliers is known in advance. 


7. Verbose outputs for debugging. 

For more details on the progress of the program and explanations of possible warnings and outcomes, set the debug parameter option to true.
