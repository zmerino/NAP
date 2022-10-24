
disp("test print new nmem")

[failed, y, pdf] = EstimatePDF(randn(1000, 1));


disp(failed)

figure('Name','Test NMEM NEW')
plot(y, pdf);

% if ~failed
%   plot(y, pdf);
% end