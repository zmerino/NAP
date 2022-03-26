function [u,sqr] = SQR(sx,sPDF,sample)

Ns = length(sample);
x = sample;

sx = sort(sx);

pL = (0.5)/Ns;          %probability for data to be  left of window
pR = (0.5)/Ns;          %probability for data to be right of window
pNorm = 1 - pL - pR;    %probability for data to fall within window

sCDF = zeros( size(sPDF) );
%----------------------------------------------------------------------------------------------------------------------------------------------
sCDF(1) = 0;
kmax = length(sCDF);
disp(['length(sCDF): ',num2str(kmax)])
%pause
for k=2:kmax
    fave = 0.5*( sPDF(k) + sPDF(k-1) );
    area = fave*( sx(k) - sx(k-1) );
    sCDF(k) = sCDF(k-1) + area;
    
    %disp(['fave: ',num2str(fave),'  area: ',num2str(area),'  sCDF(k): ',num2str(sCDF(k))])
end
temp = sCDF(kmax);
% sCDF = sCDF/temp;  % recalling what pL, pR and pNorm are
sCDF = pNorm*(sCDF/temp) + pL;  % recalling what pL, pR and pNorm are
%pL = (0.5 + nLoutliers)/Ns;    probability for data to be  left of window
%pR = (0.5 + nRoutliers)/Ns;    probability for data to be right of window
%pNorm = 1 - pL - pR;           probability for data to fall within window

sample = x(2:end-1);
%----------------------------------------------------------- adjust u range
sampleUpLim = max(sx);
sampleLoLim = min(sx);

[row, ~] = find(sample <= sampleUpLim & sample >= sampleLoLim);
% disp('sample size')
% disp(size(sample))

%u = interp1(sx,sCDF,sample);    % get corresponding u for each x in sample
u = interp1(sx,sCDF,sample(row));    % get corresponding u for each x in sample

% figure('Name','interp1 figure')
% hold on
% plot(sx,sCDF,'.r')
% plot(sample,1,'xb')
% plot(sample(row),u,'om')
% legend('sCDF','sample','u')
% pause
%--------------------------------------------------------------------------
% original
%uref = (1:Ns-2)/(Ns - 1);              % both end points have been removed
% uref = (1:Ns-2)/(Ns - 1);
uref = (1:size(sample(row),1))/(size(sample(row),1) - 1);

if( size(uref,1) ~= size(u,1) )
    u = u';
end
% --------------------------------------------------- get scaled residual

% disp(['size u: ',num2str(size(u)),' uref: ',num2str(size(uref))])
% disp(['min u: ',num2str(min(u)),' min uref: ',num2str(min(uref))])
% disp(['max u: ',num2str(max(u)),' max uref: ',num2str(max(uref))])
% disp(u)
% pause
u = sort(u);
% issue with u having NaNs: from sample range being outside of sx range
sqr = sqrt(Ns)*(u - uref); % normal formula has sqrt(Ns+2) but Ns -> Ns-2
end