x = [1,2,3,4,5,6,7,8,9]

y = [3,4,5]

% mask = x==y

bl = min(y)
br = max(y)

test1 = (min(y)<=x & max(y)>=x)

% [test3, idx] = test1(test1==1)

test4 = find(test1==1)

test1(test4(end)) = 0;

test = x(test1)

% test = find(x==y)


test3 = ismember(x,y)