% module = factor_module_iid(D, N, ...)
%
% Assume x_n ~ N(0,I)
function module = factor_module_iid(varargin)

% Parse options
options = struct('update_hyperparameters', true, ...
                 'prior', struct(), ...
                 'noise_module', [], ...
                 'init', []);
[options, errmsg] = argparse(options, varargin{:});
error(errmsg);

prior = struct();

X = [];
CovX = [];
XX = [];
rho = [];
logrho = [];

module.get_struct = @get_struct;
module.initialize = @initialize;
module.update = @update;
module.rotate = @rotate;
module.rotation_cost = @rotation_cost;

  function S = get_struct()
  % TODO: not ready yet..
  S.module = 'factor_module_iid';
  S.mu = X;
  S.Cov = CovX;
  S.rho = rho;
  S.prior = prior;
  S.options = options;
  end

  function [Z, CovZ, rhoz] = initialize(D,N)

  % Allocate memory
  X = zeros(D,N);
  CovX = zeros(D,D,N);
  rho = zeros(N,1);
  
  % Parse prior
  prior = struct('mu', 0, ...
                 'CovX', eye(D));
  [prior, errmsg] = argparse(prior, options.prior);
  error(errmsg);
  if ndims(prior.CovX)~=2 || any(size(prior.CovX)~=D)
    size(prior.CovX)
    D
    error('Prior covariance matrix no DxD matrix.');
  end
  prior.L = chol(prior.CovX, 'lower');
  prior.invCovX = linsolve_lchol(prior.L, eye(D));
  prior.logdet_CovX = logdet_chol(prior.L);
  prior.mu = prior.mu .* ones(D,1);

  % Parse custom initialization
  if ~isempty(options.init) && isstruct(options.init)
    if isfield(options.init, 'X')
      X = options.init.X;
    end
    if isfield(options.init, 'CovX')
      CovX = options.init.CovX;
    end
    if isfield(options.init, 'rho')
      rho = options.init.rho;
    end
  end

  % Default initialization
  if isempty(CovX)
    CovX = repmat(prior.CovX, [1 1 N]);
  end
  if isempty(X)
    X = mvnrnd(prior.mu, prior.CovX, N)';
  end
  if isempty(rho)
    if ~isempty(options.noise_module)
      rho = options.noise_module.initialize();
    else
      rho = ones(N,1);
    end
  end
  
  Z = X;
  CovZ = CovX;
  rhoz = rho;
  end

  function [Z, CovZ, rhoz, logrhoz, KL] = update(iter, Y, Obs, W, CovW, rhow, Tau)
  
  %
  % Update variables
  %
  
  D = size(W,1);
  [M,N] = size(Y);
  
  % Update X
  
  % Compute: Cov_n = sum_m ( tau_mn * CovW_m )
  if ndims(CovW) == 2 && all(size(CovW)==[D,M])
    Tau_CovW = CovW*Tau;
  else
    % weighted sum of covariance matrices
    M = size(W,2);
    Tau_CovW = reshape(reshape(CovW,[D*D,M])*Tau,[D,D,N]);
  end

  % Compute E[X] and Cov[X]
  V = bsxfun(@times, rhow, Tau);
  Z = zeros(size(X));
  logdet_CovX = 0;
  for n=1:N
    v = V(:,n);
    if ndims(CovW) == 2 && all(size(CovW)==[D,M])
      ww = W*spdiag(v)*W' + diag(Tau_CovW);
    else
      ww = W*spdiag(v)*W' + Tau_CovW(:,:,n);
    end
    [L,p] = chol(prior.invCovX+ww, 'lower');
    if p~=0
      v
      ww
      rcond_ww = rcond(ww)
      n
      error('Matrix not positive definite');
    end
    logdet_CovX = logdet_CovX - logdet_chol(L);
    CovX(:,:,n) = linsolve_lchol(L, eye(D));
    Z(:,n) = W*(v.*Y(:,n)) + prior.invCovX*prior.mu;
    X(:,n) = CovX(:,:,n) * Z(:,n);
  end
  
  % Update q(rho)
  if ~isempty(options.noise_module)
    %VY = V.*Y;
    E2 = dot(V.*Y,Y,1)' ...
         + prior.mu'*prior.invCovX*prior.mu ...
         - dot(X,Z,1)';
    [rho,logrho,KL_rho] = options.noise_module.update(iter, E2, sum(Obs, 1)');
  else
    rho = ones(N,1);
    logrho = zeros(N,1);
    KL_rho = 0;
  end
  
  XX = X*spdiag(rho)*X' + sum(CovX,3);
  X0X0 = XX - 2*prior.mu*(X*rho)' + sum(rho)*prior.mu*prior.mu';
  
  %
  % Compute Kullback-Leibler divergence KL(q(X,alpha)||p(X,alpha))
  %
  
  % <log p(X|alpha,rho)>
  logp_X = gaussian_logpdf(traceprod(prior.invCovX,X0X0), ...
                           0, ...
                           0, ...
                           N*prior.logdet_CovX - D*sum(logrho), ...
                           N*D);
                              
  % <log q(X|rho)>
  logq_X = -gaussian_entropy(logdet_CovX - D*sum(logrho), N*D);
  
  %KL_rho
  KL_X = -logp_X + logq_X;
  
  KL = KL_X + KL_rho;
    
  Z = X;
  CovZ = CovX;
  rhoz = rho .* ones(N,1);
  logrhoz = logrho .* ones(N,1);

  end
  
  function [c, dc] = rotation_cost(A, U, S, V)
  
  [D,N] = size(X);
  
  AX0X0A = A*XX*A' - 2*prior.mu*(X*rho)'*A' + sum(rho)*prior.mu*prior.mu';
  dAX0X0A = 2*A*XX - 2*prior.mu*(X*rho)';
  
% $$$   c = trace(AX0X0A);
% $$$   dc = dAX0X0A;
% $$$   return,
  
  invA = V*diag(1./diag(S))*U';

  c = 0;
  dc = 0;

  % Cost from <log q(X)>
  logq_X = -N * logdet_diag(S);
  dlogq_X = -N * invA'; %inv(A)';
  c = c + logq_X;
  dc = dc + dlogq_X;
  
  % Cost from -<log p(X|alpha)>
  logp_X = gaussian_logpdf(traceprod(prior.invCovX,AX0X0A), ...
                           0, ...
                           0, ...
                           N*prior.logdet_CovX - D*sum(logrho), ...
                           N*D);
  dlogp_X = gaussian_dlogpdf(prior.invCovX*dAX0X0A, ...
                             0, ...
                             0, ...
                             0);
  c = c - logp_X;
  dc = dc - dlogp_X;
  
  end
  
  function [Z, CovZ] = rotate(A)
  [D,N] = size(X);
  X = A*X;
  for n=1:N
    CovX(:,:,n) = A*CovX(:,:,n)*A';
  end
  Z = X;
  CovZ = CovX;
  end

end
