% GAUSSIAN_DLOGPDF - The log pdf of a multivariate normal distribution.
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
% X = GAUSSIAN_DLOGPDF(Y, MU, L)
%
% X = GAUSSIAN_DLOGPDF(INVCOV_Y, INVCOV_MU)
%
% X = GAUSSIAN_DLOGPDF(D_Y_INVCOV_Y, D_Y_INVCOV_MU, D_MU_INVCOV_MU, ...
%                      D_LOGDET_COV)

% Last modified 2011-02-09
% Copyright (c) Jaakko Luttinen (jaakko.luttinen@tkk.fi)

function dx = gaussian_dlogpdf(p1, p2, p3, p4)

if nargin == 3

  y = p1;
  mu = p2;
  L = p3;

  dx = -linsolve_lchol(L,y-mu);

elseif nargin == 2
  
  invCov_y = p1;
  invCov_mu = p2;

  dx = invCov_mu - invCov_y;
  
elseif nargin == 4
  
  dx = -0.5*(p1-2*p2+p3) - 0.5*p4;

end


