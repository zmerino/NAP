function val = mse(fact, fs)


if size(fact) ~= size(fs)
    error('Size of fact is different than fs')
end

n = length(fs);
val = (1/n)*sum((fact-fs).^2);

end