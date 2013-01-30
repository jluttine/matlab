% LAPLACE_LOGPDF - The log pdf of Laplace distribution.
%
%           pdf = lambda/2 * exp(-lambda*abs(x-mu))
%
% x = laplace_logpdf(x, mu, lambda)
% x = laplace_logpdf(x_lambda_mu, log_lambda)

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

function lp = laplace_logpdf(x, mu, lambda)

if nargin == 3
  x_lambda_mu = lambda * abs(x-mu);
  log_lambda = log(lambda);
else
  x_lambda_mu = x;
  log_lambda = mu;
end

lp = log_lambda - log(2) - x_lambda_mu;


