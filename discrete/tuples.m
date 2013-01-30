% Y = TUPLES(X1, X2, X3, ...)
%
% Forms the complete set of tuples from the given vectors.
%
% For instance,
%
% >> tuples([1 2], 3, [4 5 6])
%
% ans =
%
%      1     3     4
%      1     3     5
%      1     3     6
%      2     3     4
%      2     3     5
%      2     3     6
%
% The number of input vectors can be arbitrary.

% Copyright (c) 2010 Jaakko Luttinen

function Y = tuples(varargin)

% Number of lists
D = length(varargin);

N = 1;
for d=1:D
  N = N * length(varargin{d});
end

Y = zeros(N,D);

N_d = 1;
for d=D:-1:1
  l = length(varargin{d});
  ind = mod( floor((0:(N-1))/N_d), l ) + 1;
  N_d = N_d * l;
  Y(:,d) = varargin{d}(ind);
end