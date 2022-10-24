N = 1000;
data=[randn(N,1); randn(N,1)*2 + 35; randn(N,1) + 55];
data = data';

x = -10:0.1:60;
exactPDF = normpdf(x, 0, 1) + normpdf(x, 35, 2) + normpdf(x, 55, 1);
sumPdf = sum(exactPDF .* 0.1);
exactPDF = exactPDF ./ sumPdf;
d = [x', exactPDF'];    
    
PDFAnalyze(data, 'plotType', 'pdf', 'plotType', 'combined', 'plotType', 'sqr', 'distribution', d);

