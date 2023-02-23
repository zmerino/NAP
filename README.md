# Non-parametric Adaptive Parition (NAP) Desnity Estiamtion

The NAP estimator adaptively paritions a sample of univariant data, finds PDF estimates for all of the paritioned samples, and stitches the estiamtes together using a data driven method. The NAP estimator was developed for applications to large datasets with varying distributions in the data, as well as, applications to distributions that are notoriously difficult to estimate.

This method leverages a density estimator previously developed by Jenny Farmer and Donald Jacobs which is often refered to as [Non-parametric Maximum Entropy Method Density Esitmation](https://github.com/jennyfarmer/PDFAnalyze). The C++ code was compiled into MEX and published on [Mathworks](https://www.mathworks.com/matlabcentral/fileexchange/74834-multivariate-probability-density-estimation?s_tid=prof_contriblnk) for applications in MATLAB.

> **NOTE:** The NAP estimator is theoretically capable of using other methods of PDF estimation to find estiamtes for each parition. Kernal Density Estimation was previously shown to be compatible and produce equivalent results.

## Version
- version 1.0 - February 2023

## Contributors
- Zach D. Merino: zmerino2718@gmail.com
- Jenny Farmer: jfarmer@carolina.rr.com
- Donald Jacobs: djacobs1@uncc.edu

## Affiliations
- University of North Carolina at Charlotte, NC, USA
- University of Waterloo, ON, Canada



