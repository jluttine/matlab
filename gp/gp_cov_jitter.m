% GP_COV_JITTER - Adds a small constant to the diagonal of a covariance
%                 matrix in order to make it numerically positive definite.
%
% Usage:
%
%   COVFUNC_JITTERED = GP_COV_JITTER(COVFUNC)
%   COVFUNC_JITTERED = GP_COV_JITTER(COVFUNC, JITTER)
%
% COVFUNC is a covariance function to which the constant diagonal is
% added. The covariance matrix by COVFUNC must be square.
%
% The default value for JITTER is 1E-6.
%
% See also GP_COV.

% Last modified 2010-01-25
% Copyright (c) Jaakko Luttinen (jaakko.luttinen@tkk.fi)

function covfunc = gp_cov_jitter(covfunc, jitter)

[~, N1, N2] = covfunc();

if nargin < 2
  jitter = 1e-6;
end

if N1 ~= N2
  error('Jittered matrix must be square');
else
  N = N1;
end

covfunc = gp_cov_sum(covfunc, gp_cov_scale(gp_cov_delta(N), 'scale', jitter));
