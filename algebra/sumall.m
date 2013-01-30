% SUMALL - Sums all the elements of the given matrix.
%
% Usage:
%
%   Y = SUMALL(X)
%
% Evaluates Y = SUM(X(:))

% Last modified 2010-11-09
% Copyright (c) Jaakko Luttinen (jaakko.luttinen@tkk.fi)

function y = sumall(x)
y = sum(x(:));
