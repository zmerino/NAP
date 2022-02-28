function [sampleVec] = samplesVector(minSamplesExp,maxSamplesExp,dataTypeflag,step)

% to be function inputs
%--------------------------------------------------------------------------
%step = 1; %<---- can be changed to skip number of samples created
%minSamplesExp;
%maxSamplesExp;
%dataTypeflag = true; %<--- true/false integer powers of 2/real powers of 2

% Define a vector of samples to generate
%--------------------------------------------------------------------------
exponents = minSamplesExp:step:maxSamplesExp;
sampleVec = zeros(1,length(exponents));

if dataTypeflag  
    % Generates vector of samples from integer power 2
    sampleVec(1:length(exponents)) = 2.^exponents(1:length(exponents));
else
    % Generates vector of samples from real power 2
    for i = 1:length(exponents)
        n = minSamplesExp + i + rand;
        sampleVec(i) = floor(2^n);
    end
end
end