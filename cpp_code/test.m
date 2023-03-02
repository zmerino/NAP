
disp("test print new nmem")

data = randn(1000, 1)';
[failed, y, pdf] = EstimatePDF(randn(1000, 1));
% [ymv, pdfmv] = EstimatePDFmv(data);


disp(failed)

figure('Name','Test NMEM NEW')
plot(y, pdf);
% plot(ymv, pdfmv);

% if ~failed
%   plot(y, pdf);
% end