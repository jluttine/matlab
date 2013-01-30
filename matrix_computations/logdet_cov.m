% LOGDET_COV - Computes the log-determinant of a positive-definite, symmetric
% matrix SIGMA.
%
% X = LOGDET_COV(SIGMA)
%
% The function uses the Cholesky decomposition of SIGMA for stable and
% efficient computations.

% Last modified 2010-06-04
% Copyright (c) Jaakko Luttinen

function x = logdet_cov(Sigma)

x = logdet_chol(chol(Sigma));