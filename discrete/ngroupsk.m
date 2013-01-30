% X = NGROUPSK(N, K)
%
% The number of ways N*K objects can be divided into N groups of equal size
% (K objects per group).
%
% N = number of groups
% K = size of the groups
%
% X = number of different group division combinations

% Last modified 2011-01-28
% Copyright (c) Jaakko Luttinen (jaakko.luttinen@tkk.fi)

function x = groupsk(n, k)

x = factorial(n.*k) ./ ( factorial(n) .* factorial(k).^n );