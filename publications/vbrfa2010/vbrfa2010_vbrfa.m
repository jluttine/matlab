% Q = VBRFA(Y, D, ...)
%
% Variational Bayesian robust factor analysis.

% Last modified 2010-06-07
% Copyright (c) Jaakko Luttinen (jaakko.luttinen@tkk.fi)

function Q = vbrfa(Y, D, varargin)

[M,N] = size(Y);

options = struct( ...
    'init',         [],               ...
    'prior',        [],               ...
    'common_tau',   true,             ...
    'common_nu',    false,            ...
    'rotate',       1,                ...
    'update_w',     1,               ...
    'update_nu',    2,               ...
    'autosavetime', nan,              ...
    'autosavefile', 'vbrfa_autosave', ...
    'maxiter',      100,              ...
    'robustness',   'independent-t', ... % none / multivariate-t / independent-t
    'user_data',    []);

% Parse arguments
[options, errmsg] = argparse( options, varargin{:} );
error(errmsg);

% Missing values
Obs = ~isnan(Y);
[ObsM, ObsN] = find(Obs);

% Set priors
set_priors();

% Initialize
% $$$ mu = 0;   % Without this, Matlab complains about mu being a function. The
% $$$           % nested functions seem to cause this problem.
initialize();

% Fix some parameters?


% Helpful variables
Id = speye(D);
log2pi = log(2*pi);


logP = -inf;
loglikelihood = nan*zeros(options.maxiter,1);
start_time = cputime;
time = nan*zeros(options.maxiter,1);

lastsave = now;

% Do this until you permanently remove mu
mu = zeros(M,1);
v_mu = zeros(M,1);
for m=1:M
  mu(m) = 0;
  v_mu(m) = 0;
end

for iter=1:options.maxiter
  
  t = cputime;
  
  %%%%%%%%%%%%%%%%%%%%%%
  %% Update parameters
  
  %%%%%%%%%%%%%%%%%
  %
  % X BLOCK
  %
  %%%%%%%%%%%%%%%%%
  
  %% UPDATE X
  
  % q(X)
  for n=1:N
    Mmv = Obs(:,n);
    UTau = U(Mmv,n).*tau(Mmv);
    UWWTau = sum(bsxfun(@times, ...
                        reshape(U(Mmv,n),[1,1,sum(Mmv)]), ...
                        WtauW(:,:,Mmv)), ...
                 3);

    % Distribution
    CovX(:,:,n) = inv(UWWTau + Id);
    X(:,n) = CovX(:,:,n) * ( W(Mmv,:)' * ( UTau .* (Y(Mmv,n)-mu(Mmv)) ) ); 
    % Sufficient statistics
    XX(:,:,n) = CovX(:,:,n) + X(:,n)*X(:,n)';
    Sv{n} = CovX(:,:,n); % only for backward compatibility..
  end
  
  % KL-term: <log q(X)> - <log p(X)>
  KL_X = 0;
  for n=1:N
    KL_X = KL_X - gaussian_entropy(logdet_cov(CovX(:,:,n)), D);
    KL_X = KL_X - gaussian_logpdf(trace(XX(:,:,n)), ...
                                  0, ...
                                  0, ...
                                  0, ...
                                  D);
  end

  %%%%%%%%%%%%%%%%%
  %
  % NOISE BLOCK
  %
  %%%%%%%%%%%%%%%%%
  
  %% UPDATE U AND NU
  switch options.robustness

   case 'none'
    %% Normal PCA / FA
    U(Obs) = 1;
    logU(Obs) = 0;
    nu(:) = inf;
    KL_U = 0;
    
   case 'multivariate-t'
    %% Joint Student-t for the dimensions (Archambaeu et al), except X
    %% independent of U, i.e., p(X) is Gaussian.
    
    % Evaluate squared errors
    for m=1:M
      Nmv = Obs(m,:);
      
      E2(m,Nmv) = tau(m)*(Y(m,Nmv)-mu(m)).^2 ...
          + tau(m)*v_mu(m) ...
          - 2*tau(m)*(Y(m,Nmv)-mu(m)).*(W(m,:)*X(:,Nmv)) ...
          + mtraceprod(WtauW(:,:,m), XX(:,:,Nmv));
    end

    % Update nu
    if (numel(options.update_nu)==1 && iter >= options.update_nu) ...
          || any(iter==options.update_nu)

      % Common nu
      nu(:) = t_ml(sum(E2,1), nu(1), sum(Obs,1));
 %     disp('Update nu');
    end
    
    % Update U
    a_U(1,:) = (nu(1) + sum(Obs,1)) / 2;
    b_U(1,:) = (nu(1) + sum(E2,1)) / 2;
    a_U(:,:) = repmat(a_U(1,:), [M,1]);
    b_U(:,:) = repmat(b_U(1,:), [M,1]);
    U(Obs) = a_U(Obs)./b_U(Obs);
    logU(Obs) = psi(a_U(Obs)) - log(b_U(Obs));
   
    % KL-term: <log q(U)> - <log p(U)>
    KL_U = 0;
    Nmv = sum(Obs,1) > 0;
    KL_U = KL_U - sum(gamma_entropy(a_U(1,Nmv),b_U(1,Nmv)));
    KL_U = KL_U - sum(gamma_logpdf(U(1,Nmv),logU(1,Nmv),nu(1)/2,nu(1)/2));
   
   case 'independent-t'
    %% Independent Student-t for each dimension (Luttinen et al)
    
    % Evaluate squared errors (this costs A LOT!)
    for m=1:M
      Nmv = Obs(m,:);
      
      E2(m,Nmv) = tau(m)*(Y(m,Nmv)-mu(m)).^2 ...
          + tau(m)*v_mu(m) ...
          - 2*tau(m)*(Y(m,Nmv)-mu(m)).*(W(m,:)*X(:,Nmv)) ...
          + mtraceprod(WtauW(:,:,m), XX(:,:,Nmv));
    end

    % Update nu
    if (numel(options.update_nu)==1 && iter >= options.update_nu) ...
          || any(iter==options.update_nu)
      if options.common_nu
        % Common nu
        nu(:) = t_ml(E2(Obs), nu(1));
      else
        for m=1:M
          % Separate nu
          if true
            % Type II ML
            nu(m) = t_ml(E2(m,Obs(m,:)), nu(m));
          else
            % EM ML
            nv = Obs(m,:);
            func = @(lognu) (-1-log(exp(lognu)/2)+psi(exp(lognu)/2)- ...
                             mean(logU(m,nv)-U(m,nv)));
            lognu = fzero(func, log(nu(m)));
            nu(m) = exp(lognu);
          end
        end
      end
%      disp('Update nu');
    end
    
    % Update U
    for m=1:M
      Nmv = Obs(m,:);
      a_U(m,Nmv) = (nu(m)+1) / 2;
      b_U(m,Nmv) = (nu(m) + E2(m,Nmv)) / 2;
% $$$     b_U(m,Nmv) = (nu(m) + E2(m,Nmv)*tau(m)) / 2;
    end
    U(Obs) = a_U(Obs)./b_U(Obs);
    logU(Obs) = psi(a_U(Obs)) - log(b_U(Obs));
    
    % KL-term: <log q(U)> - <log p(U)>
    KL_U = 0;
    KL_U = KL_U - sum(gamma_entropy(a_U(Obs),b_U(Obs)));
    for m=1:M
      Nmv = Obs(m,:);
      KL_U = KL_U - sum(gamma_logpdf(U(m,Nmv),logU(m,Nmv),nu(m)/2,nu(m)/2));
    end
    
  end
  
  %%%%%%%%%%%%%%%%%
  %
  % W BLOCK
  %
  %%%%%%%%%%%%%%%%%
  
  %
  % Update distributions and evaluate statistics
  %
  
  % Update q(W|tau)
  for m=1:M
    Nmv = Obs(m,:);
    UXX = sum(bsxfun(@times, reshape(U(m,Nmv),[1,1,sum(Nmv)]), XX(:,:,Nmv)),3);

    CovW(:,:,m) = inv(UXX + spdiag(w));
    W(m,:) = (U(m,Nmv).*(Y(m,Nmv)-mu(m))*X(:,Nmv)') * CovW(:,:,m);
  end
  
  % Update q(tau) when q(W,tau) = q(W|tau)q(tau)
  for m=1:M
    Nmv = Obs(m,:);
    a_tau(m) = atau + 0.5 * sum(Nmv);
    b_tau(m) = btau ...
        + 0.5 * (U(m,Nmv).*Y(m,Nmv))*Y(m,Nmv)' ...
        - 0.5 * W(m,:) / CovW(:,:,m) * W(m,:)';
  end
  if options.common_tau
    % Common tau
    a_tau(:) = atau + sum(a_tau-atau);
    b_tau(:) = btau + sum(b_tau-btau);
  end
  tau = a_tau./b_tau;
  logtau = psi(a_tau) - log(b_tau);
  for m=1:M
    WtauW(:,:,m) = tau(m)*W(m,:)'*W(m,:) + CovW(:,:,m);
  end
  
  % Update w
  if (numel(options.update_w)==1 && iter >= options.update_w) || ...
        any(iter==options.update_w)
    a_w(:) = aw + 0.5*M;
    b_w(:) = bw + 0.5*diag(sum(WtauW, 3));
%    disp('Update w');
  end
  w = a_w ./ b_w;
  logw = psi(a_w) - log(b_w);
  
  %
  % Compute KL terms
  %
  
  % KL-term: <log q(W|tau)> - <log p(W|tau)>
  KL_W = 0;
  for m=1:M
    KL_W = KL_W - gaussian_entropy(-D*logtau(m) + logdet_cov(CovW(:,:,m)), D);
    KL_W = KL_W - gaussian_logpdf(diag(WtauW(:,:,m))'*w, ...
                                  0, ...
                                  0, ...
                                  -D*logtau(m) - sum(logw), ...
                                  D);
  end
  
  % KL-term: <log q(tau)> - <log p(tau)>
  KL_tau = 0;
  if options.common_tau
    KL_tau = KL_tau - gamma_entropy(a_tau(1), b_tau(1));
    KL_tau = KL_tau - gamma_logpdf(tau(1),logtau(1),atau,btau);
  else
    KL_tau = KL_tau - sum(gamma_entropy(a_tau, b_tau));
    KL_tau = KL_tau - sum(gamma_logpdf(tau,logtau,atau,btau));
  end
  
  % KL-term: <log q(w)> - <log p(w)>
  KL_w = 0;
  KL_w = KL_w - sum(gamma_entropy(a_w, b_w));
  KL_w = KL_w - sum(gamma_logpdf(w, logw, aw, bw));
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% Calculate lower bound of the log-likelihood
  
  monitor_loglikelihood = true;
  
  if monitor_loglikelihood
    
    % Cost from Y
% $$$     logpdf_Y = 0;
% $$$     for n=1:N
% $$$       for m=1:M
% $$$         if Obs(m,n)
% $$$           logpdf_Y = logpdf_Y ...
% $$$               + gaussian_logpdf(tau(m)*U(m,n)*Y(m,n)^2, ...
% $$$                                 tau(m)*U(m,n)*Y(m,n)*(W(m,:)*X(:,n)), ...
% $$$                                 U(m,n)*traceprod(WtauW(:,:,m),XX(:,:,n)), ...
% $$$                                 - logtau(m) - logU(m,n), ...
% $$$                                 1);
% $$$         end
% $$$       end
% $$$     end
% $$$     logpdf_Y

    
    % We can take advantage of the linearity of gaussian_logpdf by summing
    % before the function call:
    Yh = W*X;
    y_invCov_y = (tau(ObsM).*U(Obs))'*(Y(Obs).^2);
    y_invCov_mu = (tau(ObsM).*U(Obs))'*(Y(Obs).*Yh(Obs));
    mu_invCov_mu = 0;
    for m=1:M
      Nmv = Obs(m,:);
      mu_invCov_mu = mu_invCov_mu + ...
            traceprod(WtauW(:,:,m),wsum(XX(:,:,Nmv),U(m,Nmv),3));
% $$$       Nv = 1:N;
% $$$       for n=Nv(Obs(m,:))
% $$$         mu_invCov_mu = mu_invCov_mu + ...
% $$$             U(m,n)*traceprod(WtauW(:,:,m),XX(:,:,n));
% $$$       end
    end
    logdet_Cov = - sum(logtau(ObsM)) - sum(logU(Obs));
    dim_y = sum(Obs(:));
    
    logpdf_Y = gaussian_logpdf(y_invCov_y, ...
                               y_invCov_mu, ...
                               mu_invCov_mu, ...
                               logdet_Cov, ...
                               dim_y);
                                     
    old_logP = logP;
    
    %    logpdf_Y
    
% $$$     KL_U
% $$$     KL_X
% $$$     KL_tau
% $$$     KL_W
% $$$     KL_w
    logP = logpdf_Y - KL_U - KL_X - KL_tau - KL_W - KL_w;

% $$$     % DEBUG STUFF
% $$$     kl_Y(iter) = logpdf_Y;
% $$$     kl_U(iter) = KL_U;
% $$$     kl_X(iter) = KL_X;
% $$$     kl_W(iter) = KL_W;
% $$$     kl_tau(iter) = KL_tau;
% $$$     kl_w(iter) = KL_w;
    
    loglikelihood(iter) = logP;
  else
    loglikelihood(iter) = nan;
  end
  
  time(iter) = cputime - start_time;

  % Debugging: Check that the bound really improves
  if iter > 1 && loglikelihood(iter) < loglikelihood(iter-1)
    logP_diff = (loglikelihood(iter-1)-loglikelihood(iter))/ ...
        abs(loglikelihood(iter-1));
    warmsg = sprintf(['Lower bound decreased %.2e percents. Bug or numerical ' ...
                      'inaccuracy?'], logP_diff);
    warning(warmsg);
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% Monitoring stuff
  
  % Show progress
  fprintf('Step %d: loglikelihood=%e (%.3f seconds)\n', ...
          iter, loglikelihood(iter), cputime-t);
  
  % Check whether to save the results
  if (now-lastsave)*3600*24 >= options.autosavetime
    fprintf('Saving to %s...', options.autosavefile);
    save(options.autosavefile, 'W','CovW','w','a_w','b_w','X','CovX','Sv', ...
         'mu','v_mu','tau', 'a_tau','b_tau','U','a_U','b_U','nu', ...
         'loglikelihood', 'time', 'iter', 'options');
    lastsave = now;
    fprintf(' done.\n');
  end

  % Transformations for speeding up
  if (numel(options.rotate)==1 && iter>=options.rotate && options.rotate>=1) ...
        || any(iter==options.rotate)
    orthogonalize();
  end

end

% $$$ figure
% $$$ kl = [kl_Y(:), kl_U(:), kl_X(:), kl_W(:), kl_tau(:), kl_w(:)];
% $$$ kl = kl(6:end,:);
% $$$ kl = eminus(kl,mean(kl,1));
% $$$ plot(kl)

% Results as a struct
Q.W = W;
Q.CovW = CovW;
Q.w = w;
Q.a_w = a_w;
Q.b_w = b_w;
Q.X = X;
Q.CovX = CovX;
Q.Sv = Sv;
Q.mu = mu;
Q.v_mu = v_mu;
Q.tau = tau;
Q.a_tau = a_tau;
Q.b_tau = b_tau;
Q.U = U;
Q.a_U = a_U;
Q.b_U = b_U;
Q.nu = nu;
Q.loglikelihood = loglikelihood(1:iter);
Q.time = time(1:iter);
Q.iter = iter;
Q.options = options;


%%%%%%%%%%%%%%%%%%%%%%
%% NESTED FUNCTIONS

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function orthogonalize

% $$$   %% Move bias (approximate: mean of X to zero)
% $$$   dmu = mean(X,2);
% $$$   X = X - repmat(dmu,1,N);
% $$$   mu = mu + W*dmu;

  %% Rotate W and X

  % Whitening of X
  [Vx,D2x] = svd(sum(XX,3)/N);
  Dx = spdiag(sqrt(diag(D2x))); % convert to sparse diagonal form
  Qx = Dx \ Vx';
  Qw = full(Vx * Dx); % funny: could be sparse if size(Dx)=[1 1] ??
  W = W * Qw;
  for m=1:M
    CovW(:,:,m) = Qw' * CovW(:,:,m) * Qw;
    WtauW(:,:,m) = Qw' * WtauW(:,:,m) * Qw;
% $$$     WW(:,:,m) = W(m,:)'*W(m,:) + CovW(:,:,m);
  end

  % Orthogonalization of W 
  [Vw,Dw] = svd(sum(WtauW,3)/M);
% $$$   [Vw,Dw] = svd(sum(WW,3)/M);
  Qx = Vw' * Qx;
  Qw = Vw;
  W = W * Qw;
  for m=1:M
    CovW(:,:,m) = Qw' * CovW(:,:,m) * Qw;
    WtauW(:,:,m) = Qw' * WtauW(:,:,m) * Qw; %W(m,:)'*W(m,:) + CovW(:,:,m);
% $$$     WW(:,:,m) = W(m,:)'*W(m,:) + CovW(:,:,m);
  end
  
  % Apply rotations to X
  X = Qx * X;
  for n = 1:N
    CovX(:,:,n) = Qx * CovX(:,:,n) * Qx';
    XX(:,:,n) = X(:,n)*X(:,n)' + CovX(:,:,n);
  end
  
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%

function initialize()

  % Initialize X
  X = nan * zeros(D,N);
  Sv = cell(N,1);
  Sv(:) = {nan*zeros(D,D)};
  CovX = nan*zeros([D,D,N]); % repmat(eye(D),[1,1,N]);

  % Initialize W
  a_w = nan * ones(D,1);
  b_w = nan * ones(D,1);
  w = 1e-3 * ones(D,1);
  logw = nan * ones(D,1);
  W = orth(randn(M,D));
  CovW = repmat(eye(D),[1,1,M]);

% $$$   % Initialize mu
% $$$   v_mu = 1e0*ones(M,1);
% $$$   mu = zeros(M,1);
  
  % Initialize tau
  a_tau = nan * ones(M,1);
  b_tau = nan * ones(M,1);
  tau = 1e1 * ones(M,1);
  logtau = nan * ones(M,1);
  
  % Initialize nu and U
  nu = 1 * ones(M,1);
  a_U = nan * ones(M,N);
  b_U = nan * ones(M,N);
  U = ones(M,N);
  logU = nan * zeros(M,N);
  
  % Initialize log likelihood
  loglikelihood = -Inf;

  % Initialize mean squared errors < (y-wx-m)^2 >
  E2 = zeros(size(Y));

  % Use given initialization (size matching is checked by using ':' )
  if isstruct(options.init)
    if isfield(options.init, 'X')
      disp('Custom initialization for X.');
      warning(['Currently, X is updated first, so given initialization has ' ...
               'no effect.']);
      X(:,:) = options.init.X;
    end
    if isfield(options.init, 'CovX')
      disp('Custom initialization for CovX.');
      fprintf('Using old CovX.\n');
      warning(['Currently, X is updated first, so given initialization has ' ...
               'no effect.']);
      CovX(:,:,:) = options.init.CovX;
    end
    if isfield(options.init, 'W')
      disp('Custom initialization for W.');
      fprintf('Using old W.\n');
      W(:,:) = options.init.W;
    end
    if isfield(options.init, 'CovW')
      disp('Custom initialization for CovW.');
      fprintf('Using old CovW.\n');
      CovW(:,:,:) = options.init.CovW;
    end
    if isfield(options.init, 'w')
      disp('Custom initialization for w.');
      w(:) = options.init.w;
    end
    if isfield(options.init, 'a_w')
      disp('Custom initialization for a_w.');
      a_w(:) = options.init.a_w;
    end
    if isfield(options.init, 'b_w')
      disp('Custom initialization for b_w.');
      b_w(:) = options.init.b_w;
    end
% $$$     if isfield(options.init, 'mu')
% $$$       disp('Custom initialization for mu.');
% $$$       fprintf('Using old mu.\n');
% $$$       mu(:) = options.init.mu;
% $$$     end
% $$$     if isfield(options.init, 'v_mu')
% $$$       disp('Custom initialization for v_mu.');
% $$$       fprintf('Using old v_mu.\n');
% $$$       v_mu(:) = options.init.v_mu;
% $$$     end
    if isfield(options.init, 'tau')
      disp('Custom initialization for tau.');
      fprintf('Using old tau.\n');
      tau(:) = options.init.tau;
    end
    if isfield(options.init, 'nu')
      disp('Custom initialization for nu.');
      fprintf('Using old nu.\n');
      nu(:) = options.init.nu;
    end
    if isfield(options.init, 'U')
      disp('Custom initialization for U.');
      U(:,:) = options.init.U;
    end
  end

  % Initialize second moments for X
  XX = nan*zeros(D,D,N);
  for n=1:N
    XX(:,:,n) = X(:,n)*X(:,n)' + CovX(:,:,n);
  end

  % Initialize second moments for W
  WtauW = zeros(D,D,M);
  for m=1:M
    WtauW(:,:,m) = tau(m)*W(m,:)'*W(m,:) + CovW(:,:,m);
  end
% $$$   WW = zeros(D,D,M);
% $$$   for m=1:M
% $$$     WW(:,:,m) = W(m,:)'*W(m,:) + CovW(:,:,m);
% $$$   end
  
end

%%%%%%%%%%%%%%%%%%%%%%%%%

function set_priors()

  % Default prior
  aw = 1e-5;
  bw = 1e-5;
  atau = 1e-5;
  btau = 1e-5;
  mumu = 0;
  taumu = 1e-5;
  
  % Custom prior
  if isstruct(options.prior)
    % Prior for w
    if isfield(options.prior, 'aw')
      aw(1) = options.prior.aw;
      disp('Custom prior parameter aw');
    end
    if isfield(options.prior, 'bw')
      bw(1) = options.prior.bw;
      disp('Custom prior parameter bw');
    end
    
    % Prior for tau
    if isfield(options.prior, 'atau')
      atau(1) = options.prior.atau
      disp('Custom prior parameter atau');
    end
    if isfield(options.prior, 'btau')
      btau(1) = options.prior.btau
      disp('Custom prior parameter btau');
    end
    
    % Prior for mu
    if isfield(options.prior, 'mumu')
      mumu(1) = options.prior.mumu
      disp('Custom prior parameter mumu');
    end
    if isfield(options.prior, 'taumu')
      taumu(1) = options.prior.taumu
      disp('Custom prior parameter taumu');
    end
  end
  
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end