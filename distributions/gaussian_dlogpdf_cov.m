% GAUSSIAN_DLOGPDF_COV - The gradient of the log pdf of a multivariate
%                        normal distribution with respect to the covariance
%                        matrix.
%
% Computes the gradient of the logarithm of the multivariate normal
% probability density function N(Y|MU,COV), where Y is the random variable,
% MU is the mean and COV is the covariance. The gradient is evaluated with
% respect to the covariance matrix COV and returned in matrix form.
%
%   d(LOG(N(Y|MU,COV))) / dCOV
%
% Usage:
%
%   DX = GAUSSIAN_DLOGPDF_COV(INVCOV_Y, INVCOV_MU, INVCOV)
%
% where
%
% INVCOV_Y   : INV(COV)*Y
% INVCOV_MU  : INV(COV)*MU
% INVCOV     : INV(COV)
%
% Letting the user to compute the terms allows greater efficiency and
% flexibility.

% Last modified 2010-11-19
% Copyright (c) Jaakko Luttinen

function dx = gaussian_dlogpdf_cov(invCov_y, invCov_mu, invCov)

% $$$ if nargin == 3
% $$$   y = varargin{1};
% $$$   mu = varargin{2};
% $$$   L = varargin{3};
% $$$   invCov_y = linsolve_chol(L, y, 'lower');
% $$$   invCov_mu = linsolve_chol(L, mu, 'lower');
% $$$   invCov = linsolve_chol(L, eye(size(L)), 'lower');
% $$$ end

dx = - 0.5 * invCov ...
     + 0.5 * (invCov_y * invCov_y'  - ...
              invCov_y * invCov_mu' - ...
              invCov_mu * invCov_y' + ...
              invCov_mu * invCov_mu');


