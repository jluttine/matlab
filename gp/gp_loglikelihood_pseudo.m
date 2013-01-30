% [l, dl_f, dl_noise] = gp_loglikelihood_pseudo(y, covfunc_pseudo, theta_f,
%                                               covfunc_noise,
%                                               theta_noise)
%
% Noise covariance must be diagonal.

function [l, dl_f, dl_noise] = gp_loglikelihood_pseudo(y, covfunc_pseudo, ...
                                                  theta_f, covfunc_noise, ...
                                                  theta_noise)

% Compute the covariance matrices
if nargout < 2
  [K_pp, K_pf, k_f] = covfunc_pseudo(theta_f);
  v = covfunc_noise(theta_noise);
else
  [K_pp, K_pf, k_f, dK_pp, dK_pf, dk_f] = covfunc_pseudo(theta_f);
  [v, dv] = covfunc_noise(theta_noise);
  dV = cell(size(dv));
  for n=1:length(dV)
    if isscalar(dv{n})
      dV{n} = dv{n}*speye(length(k_f));
    elseif isvector(dv{n})
      dV{n} = spdiag(dv{n});
    else
      dV{n} = dv{n};
    end
  end
end
if isscalar(v)
  invV = @(X) X/v;
  logdetV = length(k_f) * log(v);
elseif isvector(v)
  invV = @(X) bsxfun(@times, 1./v(:), X);
  logdetV = sum(log(v));
else
% $$$   warning(['The noise covariance function returns a matrix. Make sure ' ...
% $$$            'that it is diagonal, or use a covariance function which ' ...
% $$$            'returns a scalar or a vector which is then appropriately ' ...
% $$$            'transformed to a diagonal matrix.']);
  if issparse(v)
    LD = ldlchol(v);
    invV = @(X) linsolve_ldlchol(LD, X);
    logdetV = logdet_ldlchol(LD);
  else
    L = chol(v, 'lower');
    invV = @(X) linsolve_lchol(L, X);
    logdetV = logdet_chol(L);
  end
end

% Necessary decompositions

[L_pp,p] = chol(K_pp, 'lower');
if p~=0
  disp('K_pp not positive definite, return -inf');
  theta_f
  l = -inf;
  dl_f = nan(size(theta_f));
  dl_noise = nan(size(theta_noise));
  return
end

% var(f|pseudos), assuming that V is diagonal
R = linsolve_tril(L_pp, K_pf);
var_f = k_f - dot(R,R,1)'; 

[L_Lambda,p] = chol(speye(size(K_pp)) + R*invV(R'), 'lower');
% $$$ [L_Lambda,p] = chol(K_pp + K_pf'*invV(K_pf), 'lower');
if p~=0
  clf
  theta = [theta_f(:); theta_noise(:)]
  disp('Lambda not positive definite, return -inf');
  l = -inf;
  dl_f = nan(size(theta_f));
  dl_noise = nan(size(theta_noise));
  return
end

% Compute the pseudo log-likelihood (is this only up to a constant?)
a = invV(y);
b = linsolve_tril(L_Lambda, R*a);
y_invK_y = y'*a - b'*b;
ldet = logdet_chol(L_Lambda) + logdetV;
l = gp_logpdf(y_invK_y, ldet, length(y)) ...
    - 0.5 * sum(invV(var_f));

% Compute the gradient of the pseudo loglikelihood (if requested)
if nargout >= 2
  invK_y = a - invV(R'*linsolve_triu(L_Lambda, b, true));
  %invK_y = a - invV(R'*linsolve_triu(L_Lambda, b, true));
  d = linsolve_triu(L_pp, R*invK_y, true);
  Z = linsolve_triu(L_pp, R, true);
  C_pf = invV(Z')';
  C_pp = Z*C_pf';
  W = linsolve_tril(L_Lambda, invV(R')');
  %W = linsolve_tril(L_Lambda, invV(K_fp)');
  ZW = Z*W';
  C_pf = C_pf - ZW*W;
  C_pp = C_pp - ZW*ZW';
  clear ZW;
  Z_invV = invV(Z')';
  Z_invV_Z = Z_invV*Z';
  dl_f = zeros(numel(theta_f),1);
  for n=1:length(dl_f)
    y_invK_dK_invK_y = 2*d'*dK_pf{n}*invK_y - d'*dK_pp{n}*d;
    trace_invK_dK = 2*traceprod(C_pf, dK_pf{n}) - traceprod(C_pp, dK_pp{n});
    dl_f(n) = gp_dlogpdf(y_invK_dK_invK_y, ...
                         trace_invK_dK) ...
              - 0.5 * (sum(invV(dk_f{n})) ...
                       - 2*traceprod(Z_invV, dK_pf{n}) ...
                       + traceprod(Z_invV_Z, dK_pp{n}));
  end
end

if nargout >= 3
  dl_noise = zeros(numel(theta_noise),1);
  for n=1:length(dl_noise)
    y_invK_dV_invK_y = invK_y' * dV{n} * invK_y;
    trace_invK_dV = trace(invV(dV{n})) - traceprod(W, W*dV{n});
    dl_noise(n) = gp_dlogpdf(y_invK_dV_invK_y, ...
                             trace_invK_dV) ...
        + 0.5 * sum(invV(dV{n}*invV(var_f)));
  end
end

return


% $$$ [L_Lambda,p] = chol(K_pp + K_fp'*invV(K_fp), 'lower');
% $$$ if p~=0
% $$$   clf
% $$$   theta = [theta_f(:); theta_noise(:)]
% $$$   disp('Lambda not positive definite, return -inf');
% $$$   l = -inf;
% $$$   dl_f = nan(size(theta_f));
% $$$   dl_noise = nan(size(theta_noise));
% $$$   return
% $$$ end
% $$$ 
% $$$ 
% $$$ % var(f|pseudos), assuming that V is diagonal
% $$$ Z = linsolve_tril(L_pp, K_fp');
% $$$ var_f = k_f - dot(Z,Z,1)'; 
% $$$ 
% $$$ % Compute the pseudo log-likelihood (is this only up to a constant?)
% $$$ a = invV(y);
% $$$ b = linsolve_tril(L_Lambda, K_fp'*a);
% $$$ y_invK_y = y'*a - b'*b;
% $$$ ldet = logdet_chol(L_Lambda) - logdet_chol(L_pp) + logdetV;
% $$$ l = gp_logpdf(y_invK_y, ldet, length(y)) ...
% $$$     - 0.5 * sum(invV(var_f));
% $$$ 
% $$$ % Compute the gradient of the pseudo loglikelihood (if requested)
% $$$ if nargout >= 2
% $$$   invK_y = a - invV(K_fp*linsolve_triu(L_Lambda, b, true));
% $$$   d = linsolve_lchol(L_pp, K_fp'*invK_y);
% $$$   Z = linsolve_triu(L_pp, Z, true);
% $$$   C_fp = invV(Z');
% $$$   C_pp = Z*C_fp;
% $$$   W = linsolve_tril(L_Lambda, invV(K_fp)');
% $$$   ZW = Z*W';
% $$$   C_fp = C_fp - W'*ZW';
% $$$   C_pp = C_pp - ZW*ZW';
% $$$   clear ZW;
% $$$   invV_Z = invV(Z');
% $$$   Z_invV_Z = Z*invV_Z;
% $$$   dl_f = zeros(numel(theta_f),1);
% $$$   for n=1:length(dl_f)
% $$$     y_invK_dK_invK_y = 2*invK_y'*dK_fp{n}*d - d'*dK_pp{n}*d;
% $$$     trace_invK_dK = 2*traceprod(C_fp, dK_fp{n}) - traceprod(C_pp, dK_pp{n});
% $$$     dl_f(n) = gp_dlogpdf(y_invK_dK_invK_y, ...
% $$$                          trace_invK_dK) ...
% $$$               - 0.5 * (sum(invV(dk_f{n})) ...
% $$$                        - 2*traceprod(invV_Z, dK_fp{n}) ...
% $$$                        + traceprod(Z_invV_Z, dK_pp{n}));
% $$$   end
% $$$ end
% $$$ 
% $$$ if nargout >= 3
% $$$   dl_noise = zeros(numel(theta_noise),1);
% $$$   for n=1:length(dl_noise)
% $$$     y_invK_dV_invK_y = invK_y' * dV{n} * invK_y;
% $$$     trace_invK_dV = trace(invV(dV{n})) - traceprod(W, W*dV{n});
% $$$     dl_noise(n) = gp_dlogpdf(y_invK_dV_invK_y, ...
% $$$                              trace_invK_dV) ...
% $$$         + 0.5 * sum(invV(dV{n}*invV(var_f)));
% $$$   end
% $$$ end
