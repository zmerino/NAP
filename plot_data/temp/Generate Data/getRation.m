function r = getRation(sample,window)

% disp('tic/toc start --------------------------------------------------')
% tic
sample = sort(sample);
n = length(sample);
% toc

dx = zeros(1,n-1);
% tic
dx(1:n-1) = sample(2:n) - sample(1:n-1);
% toc

% tic
dx = sort(dx');
% pause
% toc
% 
% disp('size dx: ')
% disp(size(dx))

% tic
dxMin = mean(dx(1:window));
dxMax = mean(dx(end-window:end));

% disp(['dxMin: ',num2str(dxMin)])
% disp(['dxMax: ',num2str(dxMax)])
% toc

r = dxMax/dxMin;
% disp('tic/toc end --------------------------------------------------')
end