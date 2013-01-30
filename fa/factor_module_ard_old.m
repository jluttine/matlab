% module = factor_module_ard(D, N, ...)
%
% Assume x_n ~ N(0,diag(1./alpha)) and alpha_d ~ G(a_alpha, b_alpha)

function module = factor_module_ard(D, N, varargin)

% Parse options
options = struct('update_hyperparameters', true, ...
                 'prior', struct(), ...
                 'init', []);
[options, errmsg] = argparse(options, varargin{:});
error(errmsg);

% Parse prior
prior = struct('a_alpha', 1e-6, ...
               'b_alpha', 1e-6);
[prior, errmsg] = argparse(prior, options.prior);
error(errmsg);

% The ARD variable
alpha = [];
%prior.a_alpha = prior.a_alpha.*ones(D,1);
%prior.b_alpha = prior.b_alpha.*ones(D,1);
a_alpha = nan(D,1);
b_alpha = nan(D,1);

v_X = [];
X = [];
CovX = [];
XX = [];

module.get_struct = @get_struct;
module.initialize = @initialize;
module.update = @update;
module.rotate = @rotate;
module.rotation_cost = @rotation_cost;

  function S = get_struct()
  % TODO: not ready yet..
  S.options = options;
  S.v_X = v_X;
  S.X = X;
  S.CovX = CovX;
  S.alpha = alpha;
  S.a_alpha = a_alpha;
  S.b_alpha = b_alpha;
  S.prior = prior;
  end

  function [v_Z, Z, CovZ] = initialize()
  alpha = 1e-3*ones(D,1);
  if isempty(options.init)
    X = randn(D,N);
  else
    X = options.init;
  end
  CovX = 1*ones(D,D,N);
  v_X = 1;
  
  v_Z = v_X*ones(N,1);
  Z = X;
  CovZ = CovX;
  end

  function [v_Z, Z, CovZ, KL] = update(iter, Y, Obs, v_W, W, CovW, Tau)
  
  Tau(~Obs) = 0;

  %
  % Update variables
  %
  
  % Update X
  logdet_CovX = 0;
  if ndims(CovW) == 2 && all(size(CovW)==[D,M])
    Tau_CovW = CovW*Tau;
  else
    % weighted sum of covariance matrices
    M = size(W,2);
    Tau_CovW = reshape(reshape(CovW,[D*D,M])*Tau,[D,D,N]);
  end
  for n=1:N
    taun = Tau(:,n);
    if ndims(CovW) == 2 && all(size(CovW)==[D,M])
      ww = W*spdiag(taun)*W' + diag(Tau_CovW);
    else
      ww = W*spdiag(taun)*W' + Tau_CovW(:,:,n);
    end
    [L,p] = chol(diag(alpha)+ww, 'lower');
    if p~=0
      alpha
      ww
      error('Matrix not positive definite');
    end
    logdet_CovX = logdet_CovX - logdet_chol(L);
    CovX(:,:,n) = linsolve_lchol(L, eye(D));
    X(:,n) = CovX(:,:,n) * (W*(taun.*Y(:,n)));
  end
  
  XX = X*X' + sum(CovX,3);
  xx = diag(XX);
  
  % Update alpha
  a_alpha(:) = prior.a_alpha(:) + 0.5*N;
  b_alpha(:) = prior.b_alpha(:) + 0.5*xx;
  alpha = a_alpha./b_alpha;
  logalpha = psi(a_alpha) - log(b_alpha);

  %
  % Compute Kullback-Leibler divergence KL(q(X,alpha)||p(X,alpha))
  %
  
  % <log p(X|alpha)>
  logp_X = gaussian_logpdf(alpha'*xx, ...
                           0, ...
                           0, ...
                           -N*sum(logalpha), ...
                           N*D);
                              
  % <log q(X)>
  logq_X = -gaussian_entropy(logdet_CovX, N*D);
  
  % <log p(alpha)>
  logp_alpha = sum(gamma_logpdf(alpha, prior.a_alpha, prior.b_alpha, logalpha));
  
  % <log q(alpha)>
  logq_alpha = -sum(gamma_entropy(a_alpha, b_alpha));
  
  KL = -logp_X + logq_X - logp_alpha + logq_alpha;
    
  v_Z = v_X * ones(N,1);
  Z = X;
  CovZ = CovX;

  end
  
  function [c, dc] = rotation_cost(A, U, S, V)
  %error('Not yet implemented')
  
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

  % Cost from <log q(X)>
  logq_X = -N * logdet_diag(S);
  dlogq_X = -N * invA'; %inv(A)';
  c = c + logq_X;
  dc = dc + dlogq_X;
  
  % Cost from -<log p(X|alpha)>
  logp_X = gaussian_logpdf(A_alpha'*aXXa, 0, 0, -N*sum(A_logalpha), N*D);
  dlogp_X = gaussian_dlogpdf(diag(aXXa)*dA_alpha + diag(A_alpha)*daXXa, ...
                             0, 0, -N*dA_logalpha);
  c = c - logp_X;
  dc = dc - dlogp_X;
  
  end
  
  function [Z, CovZ] = rotate(A)
  X = A*X;
  for n=1:N
    CovX(:,:,n) = A*CovX(:,:,n)*A';
  end
  XX = X*X'+sum(CovX,3);
  b_alpha = prior.b_alpha + 0.5*diag(XX);
  alpha = a_alpha ./ b_alpha;
  logalpha = psi(a_alpha) - log(b_alpha);
  Z = X;
  CovZ = CovX;
  end

end
