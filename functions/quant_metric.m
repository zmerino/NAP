function [metric_dist, quantiles, xq, fq] = quant_metric(xs, fs, Fs, qnum, d_name, metric)

actual = distributions;
actual.generate_data = false;

q_min = 0.01;
q_max = 0.99;

quantiles = linspace(q_min, q_max, qnum);
xq = interp1(Fs,xs,quantiles);
fq = interp1(xs,fs,xq);

% get actual distribution
actual.min_limit = min(xs);
actual.max_limit = max(xs);
actual.dist_name = d_name;

metric_dist = zeros(qnum,1);
for q = 2:qnum

    trunc_xs = linspace(xq(q-1), xq(q), 1000);
    trunc_fs = interp1(xs,fs,trunc_xs);

    actual.x = linspace(min(trunc_xs),max(trunc_xs),length(trunc_xs));
    actual = actual.dist_list();

    fact = actual.pdf_y;

    if isequal(upper(metric), 'MSE')
        metric_dist(q-1) = mse(fact, trunc_fs);
    elseif isequal(upper(metric), 'KL')
        metric_dist(q-1) = KLDiv(fact, trunc_fs);
    end

end

end