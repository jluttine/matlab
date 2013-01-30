% RMSEW - Evaluates weighted root mean squared error.
%
% E = RMSEW(Y,W)
%
% Y'*DIAG(W)*Y / SUM(W)

% Last modified 2011-01-11
% Copyright (c) Jaakko Luttinen (jaakko.luttinen@tkk.fi)

function e = rmsew(Y, W)

if nargin < 2
  W = ones(size(Y));
end

e = sqrt( ((W(:).*Y(:))'*Y(:)) / sum(W(:)) );
