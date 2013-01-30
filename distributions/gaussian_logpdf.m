% GAUSSIAN_LOGPDF - The log pdf of a multivariate normal distribution.
%
% Computes the logarithm of the multivariate normal probability density
% function N(Y|MU,COV), where Y is the random variable, MU is the mean
% and COV is the covariance.
%
% The function is called as:
%
% L = GAUSSIAN_LOGPDF(Y_INVCOV_Y, Y_INVCOV_MU, MU_INVCOV_MU, LOGDET_COV, D)
%
% where
%
% Y_INVCOV_Y   : Y'*INV(COV)*Y
% Y_INVCOV_MU  : Y'*INV(COV)*MU
% MU_INVCOV_MU : MU'*INV(COV)*MU
% LOGDET_COV   : LOG(DET(COV))
% D            : the dimensionality of the distribution
%
% However, these terms should be computed more efficiently than using INV
% and LOG(DET(..)).
%
% Letting the user to compute the terms allows greater efficiency and
% flexibility.
%
% Also note that the function is linear with respect to its arguments. Thus,
% some efficiency may be obtained for a large set of pdfs by summing the
% terms and calling this function only once with scalar arguments.
%
% Usage:
%
% X = GAUSSIAN_LOGPDF(Y, MU, L)
%
% X = GAUSSIAN_LOGPDF(Y_INVCOV_Y, Y_INVCOV_MU, MU_INVCOV_MU, LOGDET_COV, D)

% Last modified 2010-06-04
% Copyright (c) Jaakko Luttinen

function [x, dx] = gaussian_logpdf(y_invCov_y,y_invCov_mu,mu_invCov_mu,logdet_Cov,D)

if nargin == 3
  y = y_invCov_y;
  mu = y_invCov_mu;
  L = mu_invCov_mu;
  D = numel(y);

  z = linsolve_lchol(L,y-mu);
  x = - 0.5 * D * log(2*pi) ...
      - 0.5 * logdet_chol(L) ...
      - 0.5 * (y-mu)' * z;
  
  if nargout >= 2
    dx = -z;
  end

else
  x = - 0.5 * D * log(2*pi) ...
      - 0.5 * logdet_Cov ...
      - 0.5 * (y_invCov_y - 2*y_invCov_mu + mu_invCov_mu);
end


