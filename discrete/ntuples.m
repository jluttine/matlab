% Y = NTUPLES(X)
%
% Forms tuples for vectors 0:X(1), 0:X(2), ..., 0:X(END).
%
% For instance,
%
% >> ntuples([1 0 2])
%
% ans =
%
%      0     0     0
%      0     0     1
%      0     0     2
%      1     0     0
%      1     0     1
%      1     0     2
%
% See TUPLES for more information.

% Last modified 2011-01-28
% Copyright (c) Jaakko Luttinen (jaakko.luttinen@tkk.fi)

function Y = ntuples(x)

if ~isvector(x)
  error('Input should be a vector');
end

D = length(x);
C = cell(length(x),1);

for d=1:D
  C{d} = 0:x(d);
end

Y = tuples(C{:});

% $$$ % Total number of combinations
% $$$ N = prod(x+1);
% $$$ 
% $$$ Y = zeros(N,D);
% $$$ 
% $$$ Y(:,D) = mod(0:(N-1),x(D)+1);
% $$$ 
% $$$ for d=(D-1):-1:1
% $$$   Y(:,d) = mod( floor((0:(N-1))/prod(x((d+1):end)+1)), x(d)+1 );
% $$$ end
