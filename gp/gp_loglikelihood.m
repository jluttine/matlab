% [l, dl] = gp_loglikelihood(y, covfunc, theta)

function [l, dl] = gp_loglikelihood(y, covfunc, theta, varargin)

options = struct( ...
    'sparse2dense', true);

% Parse arguments
[options, errmsg] = argparse( options, varargin{:} );
error(errmsg);


% Compute the covariance matrix
if nargout < 2
  K = covfunc(theta);
else
  [K, dK] = covfunc(theta);
end
% $$$ if isscalar(opts.noise)
% $$$   if opts.noise ~= 0
% $$$     opts.noise = opts.noise * speye(size(K));
% $$$   end
% $$$ elseif any(size(K)~=size(opts.noise))
% $$$   error('Size mismatch between the covariance matrices K and noise.');
% $$$ end
% $$$ K = K + opts.noise;

if issparse(K) && nnz(K)/numel(K) > 0.3 && options.sparse2dense
  fprintf('Sparse to dense (density %.2f)\n', ...
          nnz(K)/numel(K));
  K = full(K);
end

% Decompose the covariance matrix
if issparse(K)
  [LD,p,q] = ldlchol(K);
  if p~=0
    error('Matrix must be positive definite.')
  end
  z = linsolve_ldlchol(LD,y,q);
  ldet = logdet_ldlchol(LD);
else
  [L,p] = chol(K, 'lower');
  if p~=0
    error('Matrix must be positive definite.')
  end
  z = linsolve_lchol(L,y);
  ldet = logdet_chol(L);
end

% Evaluate the log-likelihood
l = gp_logpdf(y'*z, ldet, length(y));

% Evaluate the gradient (if requested)
if nargout >= 2
  if issparse(K)
    invK = spinv_ldlchol(LD);
  else
    invK = inv_chol(L, 'lower');
  end
  dl = zeros(size(theta));
  for n=1:numel(dl)
    y_invK_dK_invK_y = z' * dK{n} * z;
    trace_invK_dK = traceprod(invK, dK{n});
    dl(n) = gp_dlogpdf(y_invK_dK_invK_y, trace_invK_dK);
% $$$     dl(n) = gp_dlogpdf(dK{n}, invK, z);
  end
end

% $$$ bytes = 0;
% $$$ for i=1:length(s)
% $$$   bytes = bytes + s(i).bytes;
% $$$ end
% $$$ bytes

end


