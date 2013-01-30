% module = factor_module_lds(D, N, ...)
%
% p(x_n) = N(x_n | A * x_(n-1), I)
%
% p(A) = N(a_d | 0, diag(alpha))
%
% p(alpha_d) = G(alpha_d | a_alpha, b_alpha)

function module = factor_module_lds(D, N, varargin)

%
% Parse options
%
options = struct('prior', [], ...
                 'init', [],
                 'A_module', []);
[options, errmsg] = argparse(options, varargin{:});
error(errmsg);

if isempty(options.A_module)
  A_module = factor_module_ard(D, D);
end

%
% Parse prior
%
prior = struct('a_alpha', 1e-6, ...
               'b_alpha', 1e-6, ...
               'mu_x0',   zeros(D,1), ...
               'Cov_x0',  1e6*eye(D));
[prior, errmsg] = argparse(prior, options.prior);
error(errmsg);

% Innovation covariance matrix
% Q = speye(D);

%
% Parse initialization of the posterior
%
init = struct('a_alpha', [], ...
              'b_alpha', [], ...
              'mu_A',    [], ...
              'Cov_A',   [], ...
              'mu_x0',   [], ...
              'Cov_x0',  [], ...
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

% A
if ~isempty(init.mu_A)
  mu_A = broadcast(init.mu_A, [D D]);
else
  mu_A = diag(alpha.^(-0.5)) * randn(D,D);
end
if ~isempty(init.Cov_A)
  Cov_A = broadcast(init.Cov_A, [D D D]);
else
  Cov_A = repmat(diag(1./alpha), [1 1 D]);
end

% x0
if ~isempty(init.mu_x0)
  mu_x0 = broadcast(init.mu_x0(:), [D 1]);
else
  mu_x0 = broadcast(prior.mu_x0(:), [D 1]);
end
if ~isempty(init.Cov_x0)
  Cov_x0 = broadcast(init.Cov_x0, [D D]);
else
  Cov_x0 = broadcast(prior.Cov_x0, [D D]);
end

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

% Ignore rho
rho = ones(N,1);
logrho = zeros(N,1);

init = struct('a_alpha', a_alpha, ...
              'b_alpha', b_alpha, ...
              'mu_A',    mu_A, ...
              'Cov_A',   Cov_A, ...
              'mu_x0',   mu_x0, ...
              'Cov_x0',  Cov_x0, ...
              'mu_X',    mu_X, ...
              'Cov_X',   Cov_X);


module.get_struct = @get_struct;
module.update = @update;
module.rotate = @rotate;
module.rotation_cost = @rotation_cost;

  function S = get_struct()
  S.module = 'factor_module_lds';
  
  % Posterior distribution parameters
  S.posterior = struct('a_alpha', a_alpha, ...
                       'b_alpha', b_alpha, ...
                       'mu_A',    mu_A, ...
                       'Cov_A',   Cov_A, ...
                       'mu_x0',   mu_x0, ...
                       'Cov_x0',  Cov_x0, ...
                       'mu_X',    mu_X, ...
                       'Cov_X',   Cov_X);

  % Posterior expectations
  S.alpha = a_alpha./b_alpha;
  S.A = mu_A;
  S.x0 = mu_x0;
  S.X = mu_X;
  S.CovX = Cov_X;
  S.rho = rho;

  S.prior = prior;
  S.init = init;
  S.options = options;
  end

  function [Z, CovZ, rhoz, logrhoz, KL] = update(iter, Y, Obs, W, CovW, ...
                                                 rhow, Tau)

  
  M = size(Y,1);
  
  %
  % Update variables
  %
  
  % Compute: Cov_n = sum_m ( tau_mn * CovW_m )
  if ndims(CovW) == 2 && all(size(CovW)==[D,M])
    Tau_CovW = CovW*Tau;
  else
    % weighted sum of covariance matrices
    M = size(W,2);
    Tau_CovW = reshape(reshape(CovW,[D*D,M])*Tau,[D,D,N]);
  end
  
  % Update q(X,x0)

  % Kalman filtering for the augmented model
  %S_A = sum(Cov_A,3);
  U_A = sqrt_cov(sum(Cov_A,3))';
  mu_xn = prior.mu_x0;
  Cov_xn = prior.Cov_x0;
  for n=1:N
    obs = Obs(:,n);
    % Augmented Y
    yn = [Y(obs,n); zeros(2*D,1)];
    % Augmented W
    if ndims(CovW) == 2 && all(size(CovW)==[D,M])
      U_W = diag(sqrt(Tau_CovW(:,n)))';
    else
      U_W = sqrt_cov(Tau_CovW(:,:,n))';
    end
    if n < N
      Wn = [W(:,obs)'; U_W; U_A];
    else
      Wn = [W(:,obs)'; U_W; zeros(D)];
    end
    % Augmented Tau, (inv(chol(R,'lower'))
    taun = diag(sqrt([Tau(obs,n); ones(2*D,1)]));
    % Kalman filter step
    [mu_xn, Cov_xn] = kalman_filter_step(mu_xn, ...
                                         Cov_xn, ...
                                         yn, ...
                                         mu_A, ...
                                         eye(D), ...
                                         Wn, ...
                                         taun, ...
                                         2);
    % Store results
    mu_X(:,n) = mu_xn;
    Cov_X(:,:,n) = Cov_xn;
    
  end
  
  % RTS smoothing
  [mu_X, Cov_X, mu_x0, Cov_x0, Cov_X_X, entropy_X] = ...
      rts_smoother(prior.mu_x0, ...
                   prior.Cov_x0, ...
                   mu_X, ...
                   Cov_X, ...
                   mu_A, ...
                   eye(D));
  
  [mu_A, Cov_A, rho_A, logrho_A, KL_A] = ...
      update(iter, ...
             Y, ...
             Obs, ...
             W, ...
             CovW, ...
             rhow, ...
             Tau);
  
% $$$   % Compute E[X] and Cov[X]
% $$$   V = bsxfun(@times, rhow, Tau);
% $$$   for n=1:N
% $$$     v = V(:,n);
% $$$     if ndims(CovW) == 2 && all(size(CovW)==[D,M])
% $$$       ww = W*spdiag(v)*W' + diag(Tau_CovW(:,n));
% $$$     else
% $$$       ww = W*spdiag(v)*W' + Tau_CovW(:,:,n);
% $$$     end
% $$$     [L,p] = chol(diag(alpha)+ww, 'lower');
% $$$     if p~=0
% $$$       %rhow1 = rhow(1)
% $$$       alpha
% $$$       %ww
% $$$       error('Matrix not positive definite');
% $$$     end
% $$$     logdet_CovX = logdet_CovX - logdet_chol(L);
% $$$     CovX(:,:,n) = linsolve_lchol(L, eye(D));
% $$$     X(:,n) = CovX(:,:,n) * (W*(v.*Y(:,n)));
% $$$   end
% $$$   
% $$$   % Update q(rho)
% $$$   rho = ones(N,1);
% $$$   logrho = zeros(N,1);
% $$$   KL_rho = 0;
  
  
% $$$   XX = X*spdiag(rho)*X' + sum(CovX,3);
% $$$   xx = diag(XX);
  
% $$$   % Update alpha
% $$$   if index_selected(iter, options.update_alpha)
% $$$     a_alpha(:) = prior.a_alpha(:) + 0.5*N;
% $$$     b_alpha(:) = prior.b_alpha(:) + 0.5*xx;
% $$$     alpha = a_alpha./b_alpha;
% $$$     logalpha = psi(a_alpha) - log(b_alpha);
% $$$   end
  
  AA = mu_A'*mu_A + sum(Cov_A,3);
  
  %
  % Compute Kullback-Leibler divergence KL(q(X,alpha)||p(X,alpha))
  %
  
  % <log p(X|A)>
  xm = mu_x0;
  xmxm = xm*xm' + Cov_x0;
  logp_X = gaussian_logpdf(traceprod(xmxm, inv(prior.Cov_x0)), ...
                           xm' * (prior.Cov_x0 \ prior.mu_x0), ...
                           prior.mu_x0' * (prior.Cov_x0 \ prior.mu_x0), ...
                           logdet_cov(prior.Cov_x0), ...
                           D);
  for n=1:N
    xnxn = mu_X(:,n)*mu_X(:,n)' + Cov_X(:,:,n);
    xmxn = xm*mu_X(:,n)' + Cov_X_X(:,:,n);
    logp_X = logp_X + ...
             gaussian_logpdf(traceprod(xnxn, eye(D)), ...
                             traceprod(xmxn, mu_A'), ...
                             traceprod(xmxm, AA), ...
                             0, ...
                             D);
    xm = mu_X(:,n);
    xmxm = xnxn;
  end
                              
  % <log q(X)>
  logq_X = -entropy_X;
  
% $$$   % <log p(alpha)>
% $$$   logp_alpha = sum(gamma_logpdf(alpha, prior.a_alpha, prior.b_alpha, logalpha));
% $$$   
% $$$   % <log q(alpha)>
% $$$   logq_alpha = -sum(gamma_entropy(a_alpha, b_alpha));
  
% $$$   KL = -logp_X + logq_X - logp_alpha + logq_alpha + KL_rho;

  KL = -logp_X + logq_X;
    
  Z = mu_X;
  CovZ = Cov_X;
  rhoz = rho .* ones(N,1);
  logrhoz = logrho .* ones(N,1);

  end
  
  %function kf_step()
  
  function [c, dc] = rotation_cost(A, U, S, V)
  %error('Not yet implemented')
  
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
