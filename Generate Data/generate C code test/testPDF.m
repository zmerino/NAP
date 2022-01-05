clear all;

sPDF = importdata(['test.dat']);
sx = importdata(['testX.dat']);
% 
% [row, col] = find(~isnan(sPDF));
% 
% sPDF = sPDF(row);
% sx = sx(row);

pNorm = 1;
pL = 0;

sCDF = zeros( size(sPDF) );
%----------------------------------------------------------------------
sCDF(1) = 0;
kmax = length(sCDF);
% disp(['length(sCDF): ',num2str(kmax)])
for k=2:kmax
    fave = 0.5*( sPDF(k) + sPDF(k-1) );
    area = fave*( sx(k) - sx(k-1) );
    sCDF(k) = sCDF(k-1) + area;
%     disp(['sCDF: ',num2str(sCDF(k))])
%     disp(['sPDF: ',num2str(sPDF(k))])
%     if sCDF(k) ~= 0
%         pause
%     end
end
temp = sCDF(kmax);
sCDF = pNorm*(sCDF/temp) + pL;  % recalling what pL, pR and pNorm are
