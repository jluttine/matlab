% X = NPAIRSK(N, K)
%
% Number of pair combinations for N objects to K pairs.
%
% For instance,
%
%   npairsk(2,1) = 1
% 
% because two elements can be divided into a pair in one way.
%
%   npairsk(4,2) = 3
%
% because four elements (A,B,C,D) can be divided into two pairs in three
% ways: [(A,B);(C,D)], [(A,C);(B,D)], [(A,D);(B,C)].
%
%   npairsk(3,1) = 3
%
% because [(A,B)], [(A,C)], [(B,C)].

% Last modified 2011-01-28
% Copyright (c) Jaakko Luttinen (jaakko.luttinen@tkk.fi)

function x = npairsk(n, k)

x = factorial(n) ./ (factorial(2*k) .* factorial(n-2*k)) .* ngroupsk(k,2);