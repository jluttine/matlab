% module = factor_module_ard(D, N, ...)
%
% Assume x_n ~ N(0,diag(1./alpha)) and alpha_d ~ G(a_alpha, b_alpha)

function module = factor_module_ard(D, N, varargin)

%
% Parse options
%
options = struct('update_alpha', true, ...
                 'prior', struct(), ...
                 'noise_module', [], ...
                 'init', []);
[options, errmsg] = argparse(options, varargin{:});
error(errmsg);

%
% Parse prior
%
prior = struct('a_alpha', 1e-6, ...
               'b_alpha', 1e-6);
[prior, errmsg] = argparse(prior, options.prior);
error(errmsg);

% Allocate memory
% $$$ XX = [];
% $$$ logrho = [];
% $$$ a_alpha = nan(D,1);
% $$$ b_alpha = nan(D,1);

%
% Parse initialization
%
init = struct('a_alpha', [], ...
              'b_alpha', [], ...
              'mu_X',    [], ...
              'Cov_X',   []);
[init, errmsg] = argparse(init, options.init);
error(errmsg);

% alpha
if ~isempty(init.a_alpha)
  a_alpha = broadcast(init.a_alpha(:), [D 1]);
else
  a_alpha = broadcast(prior.a_alpha(:), [D 1]);
end
if ~isempty(init.b_alpha)
  b_alpha = broadcast(init.b_alpha(:), [D 1]);
else
  b_alpha = broadcast(prior.b_alpha(:), [D 1]);
end
alpha = a_alpha./b_alpha;
logalpha = psi(a_alpha) - log(b_alpha);

% X
if ~isempty(init.mu_X)
  mu_X = broadcast(init.mu_X, [D N]);
else
  mu_X = randn(D,N);
end
if ~isempty(init.Cov_X)
  Cov_X = broadcast(init.Cov_X, [D D N]);
else
  Cov_X = repmat(eye(D), [1 1 N]);
end

% rho
if ~isempty(options.noise_module)
  rho_struct = options.noise_module.get_struct();
  rho = rho_struct.Tau;
else
  rho = ones(N,1);
  logrho = zeros(N,1);
end

% $$$ %
% $$$ % Initialize
% $$$ %
% $$$ if isfield(options.init, 'X')
% $$$   X = options.init.X;
% $$$ else
% $$$   X = randn(D,N);
% $$$ end
% $$$ if isfield(options.init, 'CovX')
% $$$   CovX = options.init.CovX;
% $$$ else
% $$$   CovX = repmat(eye(D), [1 1 N]);
% $$$ end
% $$$ if isfield(options.init, 'a_alpha')
% $$$   a_alpha(:) = options.init.a_alpha;
% $$$ else
% $$$   a_alpha(:) = prior.a_alpha;
% $$$ end
% $$$ if isfield(options.init, 'b_alpha')
% $$$   b_alpha(:) = options.init.b_alpha;
% $$$ else
% $$$   b_alpha(:) = prior.b_alpha;
% $$$ end
% $$$ alpha = a_alpha./b_alpha;
% $$$ logalpha = psi(a_alpha) - log(b_alpha);
% $$$ if ~isempty(options.noise_module)
% $$$   rho_struct = options.noise_module.get_struct();
% $$$   rho = rho_struct.Tau;
% $$$ else
% $$$   rho = ones(N,1);
% $$$ end
% $$$ init = struct('X', X, ...
% $$$               'CovX', CovX, ...
% $$$               'alpha', alpha);
init = struct('a_alpha', a_alpha, ...
              'b_alpha', b_alpha, ...
              'mu_X',    mu_X, ...
              'Cov_X',   Cov_X);

module.get_struct = @get_struct;
module.update = @update;
module.rotate = @rotate;
module.rotation_cost = @rotation_cost;

  function S = get_struct()
  S.module = 'factor_module_ard';
  S.posterior = struct('a_alpha', a_alpha, ...
                       'b_alpha', b_alpha, ...
                       'mu_X',    mu_X, ...
                       'Cov_X',   Cov_X);
  S.X = mu_X;
  S.CovX = Cov_X;
  S.rho = rho;
  S.alpha = alpha;
  if ~isempty(options.noise_module)
    S.rho_struct = options.noise_module.get_struct(); 
  end
  S.prior = prior;
  S.init = init;
  S.options = options;
  end

  function [Z, CovZ, rhoz, logrhoz, KL] = update(iter, Y, Obs, W, CovW, ...
                                                 rhow, Tau)

  
  M = size(Y,1);
  %D = size(W,1);
  
  %
  % Update variables
  %
  
  % Update q(X|rho)
  
  % Compute: Cov_n = sum_m ( tau_mn * CovW_m )
  logdet_CovX = 0;
% $$$   size(Tau)
% $$$   size(Y)
% $$$   size(CovW)
  if ndims(CovW) == 2 && all(size(CovW)==[D,M])
    Tau_CovW = CovW*Tau;
  else
    % weighted sum of covariance matrices
    M = size(W,2);
    Tau_CovW = reshape(reshape(CovW,[D*D,M])*Tau,[D,D,N]);
  end
  
  % Compute E[X] and Cov[X]
  V = bsxfun(@times, rhow, Tau);
  for n=1:N
    v = V(:,n);
    if ndims(CovW) == 2 && all(size(CovW)==[D,M])
      ww = W*spdiag(v)*W' + diag(Tau_CovW(:,n));
    else
      ww = W*spdiag(v)*W' + Tau_CovW(:,:,n);
    end
    [L,p] = chol(diag(alpha)+ww, 'lower');
    if p~=0
      %rhow1 = rhow(1)
      alpha
      %ww
      error('Matrix not positive definite');
    end
    logdet_CovX = logdet_CovX - logdet_chol(L);
    Cov_X(:,:,n) = linsolve_lchol(L, eye(D));
    mu_X(:,n) = Cov_X(:,:,n) * (W*(v.*Y(:,n)));
  end
  CovX = Cov_X;
  X = mu_X;
  
  % Update q(rho)
  if ~isempty(options.noise_module)
    VY = V.*Y;
    E2 = dot(VY,Y,1)' - dot(W'*X,VY,1)'; % ??
    [rho,logrho,KL_rho] = options.noise_module.update(iter, E2, sum(Obs, 1)');
  else
    rho = ones(N,1);
    logrho = zeros(N,1);
    KL_rho = 0;
  end
  
  
  XX = X*spdiag(rho)*X' + sum(CovX,3);
  xx = diag(XX);
  
  % Update alpha
  if index_selected(iter, options.update_alpha)
    a_alpha(:) = prior.a_alpha(:) + 0.5*N;
    b_alpha(:) = prior.b_alpha(:) + 0.5*xx;
    alpha = a_alpha./b_alpha;
    logalpha = psi(a_alpha) - log(b_alpha);
  end
  
  %
  % Compute Kullback-Leibler divergence KL(q(X,alpha)||p(X,alpha))
  %
  
  % <log p(X|alpha,rho)>
  logp_X = gaussian_logpdf(alpha'*xx, ...
                           0, ...
                           0, ...
                           -N*sum(logalpha) - D*sum(logrho), ...
                           N*D);
                              
  % <log q(X|rho)>
  logq_X = -gaussian_entropy(logdet_CovX - D*sum(logrho), N*D);
  
  % <log p(alpha)>
  logp_alpha = sum(gamma_logpdf(alpha, prior.a_alpha, prior.b_alpha, logalpha));
  
  % <log q(alpha)>
  logq_alpha = -sum(gamma_entropy(a_alpha, b_alpha));
  
  KL = -logp_X + logq_X - logp_alpha + logq_alpha + KL_rho;
    
  Z = X;
  CovZ = CovX;
  rhoz = rho .* ones(N,1);
  logrhoz = logrho .* ones(N,1);

  end
  
  function [c, dc] = rotation_cost(A, U, S, V)
  %error('Not yet implemented')
  
  CovX = Cov_X;
  X = mu_X;
  XX = X*spdiag(rho)*X' + sum(CovX,3);

  [D,N] = size(X);

% $$$   %AX = A*X;
% $$$   CovAX = 0;
% $$$   for n=1:N
% $$$     CovAX = CovAX + A*CovX(:,:,n)*A';
% $$$   end
% $$$   % <rho*X*X'>
% $$$   AXXA = A*(X*spdiag(rho)*X')*A' + CovAX;

  AXXA = A*XX*A';
  aXXa = diag(AXXA);
  daXXa = 2*A*XX;
  invA = V*diag(1./diag(S))*U';

  c = 0;
  dc = 0;
  
  A_b_alpha = prior.b_alpha + 0.5*aXXa;
  A_logb_alpha = log(A_b_alpha);
  A_alpha = a_alpha ./ A_b_alpha;
  A_logalpha = psi(a_alpha) - A_logb_alpha;
  dA_b_alpha = 0.5*daXXa;
  dA_logb_alpha = diag(1./A_b_alpha) * dA_b_alpha;
  dA_alpha = diag(-a_alpha.*A_b_alpha.^(-2)) * dA_b_alpha;
  dA_logalpha = -dA_logb_alpha; % diag(-1./A_b_alpha) * dA_b_alpha;
  
  % Cost from <log q(alpha)>
  logq_alpha = -sum(gamma_entropy(a_alpha, A_b_alpha));
  dlogq_alpha = -(-dA_logb_alpha);
  c = c + logq_alpha;
  dc = dc + dlogq_alpha;

  % Cost from -<log p(alpha)>
  logp_alpha = sum(gamma_logpdf(A_alpha, prior.a_alpha, prior.b_alpha, ...
                                A_logalpha));
  dlogp_alpha = gamma_dlogpdf(dA_alpha, prior.a_alpha, prior.b_alpha, ...
                              dA_logalpha);
  c = c - logp_alpha;
  dc = dc - dlogp_alpha;

  % Cost from <log q(X|rho)>
  logq_X = -N * logdet_diag(S);
  dlogq_X = -N * invA'; %inv(A)';
  c = c + logq_X;
  dc = dc + dlogq_X;
  
  % Cost from -<log p(X|alpha)>
  logp_X = gaussian_logpdf(A_alpha'*aXXa, ...
                           0, ...
                           0, ...
                           -N*sum(A_logalpha) - D*sum(logrho), ...
                           N*D);
  dlogp_X = gaussian_dlogpdf(diag(aXXa)*dA_alpha + diag(A_alpha)*daXXa, ...
                             0, 0, -N*dA_logalpha);
  c = c - logp_X;
  dc = dc - dlogp_X;
  
  end
  
  function [Z, CovZ] = rotate(A)
  CovX = Cov_X;
  X = mu_X;

  [D,N] = size(X);
  X = A*X;
  for n=1:N
    CovX(:,:,n) = A*CovX(:,:,n)*A';
  end
  XX = X*spdiag(rho)*X'+sum(CovX,3);
  b_alpha = prior.b_alpha + 0.5*diag(XX);
  alpha = a_alpha ./ b_alpha;
  logalpha = psi(a_alpha) - log(b_alpha);
  Z = X;
  CovZ = CovX;
  end

end
