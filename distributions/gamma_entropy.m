% GAMMA_ENTROPY - Computes the entropy of a Gamma distribution.
%
% H = GAMMA_ENTROPY(A, B)
%
% where A is the shape and B is the inverse scale (note the difference to
% MATLAB's built-in functions).

% Last modified 2010-06-07
% Copyright (c) Jaakko Luttinen

function H = gamma_entropy(a,b)

H = a - log(b) + gammaln(a) - (a-1).*psi(a);

