% GAUSSIAN_ENTROPY - Computes the entropy of a multivariate normal
% distribution with covariance SIGMA.
%
% H = GAUSSIAN_ENTROPY(Y, D)
%
% Y is the log-determinant of the covariance matrix, i.e., LOG(DET(SIGMA)).
% However, it should be computed more efficiently using, e.g., LOGDET_COV.
%
% D is the dimensionality of the distribution.
%
% Also note that the function is linear with respect to its arguments. Thus,
% some efficiency may be obtained for a large set of distributions by
% summing the terms and calling this function only once with scalar
% arguments.
%
% See also LOGDET_COV.

% Last modified 2010-06-09
% Copyright (c) Jaakko Luttinen

function H = gaussian_entropy(logdet_Cov, D)

H = 0.5 * (D*log(2*pi) + logdet_Cov + D);