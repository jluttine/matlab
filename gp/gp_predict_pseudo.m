% GP_PREDICT_PSEUDO - Computes the predictive distribution of a Gaussian
%                     process using pseudo inputs.
%
% Usage:
%
%   F_MEAN = GP_PREDICT_PSEUDO(Y, V, K_FP, K_PP, K_PH, K_H)
%   [F_MEAN, F_VAR] = GP_PREDICT_PSEUDO(...)
%
% The marginal of a conditional model
%
%   Y|F ~ N(F, V)
%
% is approximated by
%
%   Y ~ N(0, K_FP*INV(K_PP)*K_FP' + V)
%
% Y : Nx1 vector of observations
% V : NxN diagonal matrix of noise covariance, or a function V(X) which
%     solves V\X
% K_FP : NxM covariance matrix of function values at observation inputs
%        and at pseudo inputs 
% K_PP : MxM covariance matrix of function values at pseudo inputs
% K_PH : MxK covariance matrix of function values at pseudo inputs and at
%        inputs to be predicted
% K_H : Kx1 vector of variances of function values at inputs to be
%       predicted
%
% See also GP_LEARN_PSEUDO, GP_PREDICT.

% Last modified 2011-01-27
% Copyright (c) Jaakko Luttinen (jaakko.luttinen@tkk.fi)

function [f_mean, f_var] = gp_predict_pseudo(y, V, K_pf, K_pp, K_ph, k_h)

if isnumeric(V)
  % Get a solver for V\z
  inv_V = get_linsolve_cov(V);
else
  % A function handle which solves V\z
  inv_V = V;
end

% TODO: 
%
% - Take K_pf instead of K_fp
%
% - Take V as a vector?

[L_p,p] = chol(K_pp, 'lower');
if p~=0
  figure
  imagesc(K_pp);
  error('Matrix not positive definite');
end

Z_f = linsolve_tril(L_p, K_pf);

Lambda = speye(size(K_pp)) + Z_f*inv_V(Z_f');
inv_Lambda = get_linsolve_cov(Lambda);

f_mean = K_ph' * linsolve_triu(L_p, inv_Lambda(Z_f*inv_V(y)), true);
if nargout >= 2
  Z_h = linsolve_tril(L_p, K_ph);
  f_var = k_h(:) - dot(Z_h, Z_h - inv_Lambda(Z_h),1)';
end

% $$$ Lambda = K_pp + K_fp' * inv_V(K_fp);
% $$$ inv_Lambda = get_linsolve_cov(Lambda);
% $$$ inv_Kpp = get_linsolve_cov(K_pp);
% $$$ f_mean = K_ph' * inv_Lambda(K_fp' * inv_V(y));
% $$$ f_var = k_h(:) - dot(K_ph, inv_Kpp(K_ph) - inv_Lambda(K_ph),1)';
