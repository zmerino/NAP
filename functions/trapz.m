function Fs = trapz(xs,fs)

% length of sample estimate
ns = length(xs);
% numerically integrate PDF
Fs = zeros(ns,1);
for n = 2:ns
    Fs(n) = Fs(n-1)+ (xs(n)-xs(n-1))*(fs(n) + fs(n-1))/2;
end

end