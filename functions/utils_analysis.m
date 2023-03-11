classdef utils_analysis
    methods(Static)
        function [htx, hty, htx_act, hty_act] = heavy(xcdf, cdf, obj)

            basex = 1;
            %     basey = log(1000);
            basey = 1;

            ind = find(xcdf>0);

            htx = log(xcdf(ind)) / basex;
            hty = log(1-cdf(ind)) / basey;

            obj.x = xcdf(ind);
            obj = dist_list(obj);

            htx_act = log(obj.x) / basex;
            hty_act = log(1-obj.cdf_y) / basey;
        end

        function dist=kl(P,Q)
            %  dist = KLDiv(P,Q) Kullback-Leibler divergence of two discrete probability
            %  distributions
            %  P and Q  are automatically normalised to have the sum of one on rows
            % have the length of one at each
            % P =  n x nbins
            % Q =  1 x nbins or n x nbins(one to one)
            % dist = n x 1

            if size(P,2)~=size(Q,2)
                if (size(P,2)==size(Q,1) && size(P,1)==size(Q,2)) || (size(P,1)==size(Q,2) && size(P,2)==size(Q,1))
                    P = P';
                else
%                     disp('P')
%                     size(P)
%                     disp('Q')
%                     size(Q)
                    error('the number of columns in P and Q should be the same');
                end
            end
            if sum(~isfinite(P(:))) + sum(~isfinite(Q(:)))
%                 sum(~isfinite(P(:)))
%                 sum(~isfinite(Q(:)))
%                 fail1 = P';
%                 fail2 = Q';
%                 fail1(~isfinite(fail1))
%                 fail2(~isfinite(fail2))
%                 min(P)
%                 max(P)
%                 min(Q)
%                 max(Q)
                error('the inputs contain non-finite values!')
            end
            % normalizing the P and Q
            if size(Q,1)==1
                t6 = Q;
                t7 = P;
                t8 = sum(Q);
                t9 = repmat(sum(P,2),[1 size(P,2)]);

                Q = Q ./sum(Q);
                P = P ./repmat(sum(P,2),[1 size(P,2)]);
                temp =  P.*log(P./repmat(Q,[size(P,1) 1]));
                temp(isnan(temp))=0;% resolving the case when P(i)==0
                dist = sum(temp,2);

                if ~isreal(dist)
                    disp('KLDiv:')
                    t1 = Q
                    t2 = P
                    t3 = log(P./repmat(Q,[size(P,1) 1]))
                    t4 = P./repmat(Q,[size(P,1) 1])
                    t5 = repmat(Q,[size(P,1) 1])
                    t6 = t6
                    t7 = t7
                    t8 = t8
                    t9 = t9
                    pause
                end
            elseif size(Q,1)==size(P,1)

                Q = Q ./repmat(sum(Q,2),[1 size(Q,2)]);
                P = P ./repmat(sum(P,2),[1 size(P,2)]);
                temp =  P.*log(P./Q);
                temp(isnan(temp))=0; % resolving the case when P(i)==0
                dist = sum(temp,2);
                %
                %     disp(['Q: ', num2str(min(repmat(Q,[size(P,1) 1]))), ' ', num2str(max(repmat(Q,[size(P,1) 1])))])
                %     disp('size Q = size P')
                %     sum(~isfinite(Q ))
                %     sum(~isfinite(P))
                %     sum(~isfinite(temp))
                %     sum(~isfinite(dist))

                if ~isreal(dist)
                    t1 = Q
                    t2 = P
                    t3 = log(P./repmat(Q,[size(P,1) 1]))
                    t4 = P./repmat(Q,[size(P,1) 1])
                    t5 = repmat(Q,[size(P,1) 1])
                    t6 = t6
                    t7 = t7
                    t8 = t8
                    t9 = t9
                    pause
                end
            end
        end

        function val = mse(fact, fs)


            if size(fact) ~= size(fs)
                error('Size of fact is different than fs')
            end

            n = length(fs);
            val = (1/n)*sum((fact-fs).^2);

        end

        function [metric_dist, quantiles, xq, fq] = quant_metric(xs, fs, Fs, qnum, d_name, metric)

            actual = distributions;
            actual.generate_data = false;

            q_min = 0.01;
            q_max = 0.99;

            quantiles = linspace(q_min, q_max, qnum);

            [Fs, a1, b1] = unique(Fs);
            xs = xs(a1);
            fs = fs(a1);

            % try
            %     xq = interp1(Fs,xs,quantiles);
            % catch
            %
            %     [test1, a1, b1] = unique(xs);
            %     [test2, a2, b2] = unique(Fs);
            %     test = 'test'
            % end

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
                    metric_dist(q-1) = utils_analysis.mse(fact, trunc_fs);
                elseif isequal(upper(metric), 'KL')
                    metric_dist(q-1) = utils_analysis.kl(fact, trunc_fs);
                end

            end

        end

        function Fs = trapz(xs,fs)

            % length of sample estimate
            ns = length(xs);
            % numerically integrate PDF
            Fs = zeros(ns,1);
            for n = 2:ns
                Fs(n) = Fs(n-1)+ (xs(n)-xs(n-1))*(fs(n) + fs(n-1))/2;
            end

        end

    end
end