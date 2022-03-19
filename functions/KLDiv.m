function dist=KLDiv(P,Q)
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
        disp('P')
        size(P)
        disp('Q')
        size(Q)
        error('the number of columns in P and Q should be the same');
    end
end
if sum(~isfinite(P(:))) + sum(~isfinite(Q(:)))
    sum(~isfinite(P(:)))
    sum(~isfinite(Q(:)))
    fail1 = P';
    fail2 = Q';
    fail1(~isfinite(fail1))
    fail2(~isfinite(fail2))
    min(P)
    max(P)
    min(Q)
    max(Q)
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






