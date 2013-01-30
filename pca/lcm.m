% LCM - Variational Bayesian inference for latent component models.

% Likelihood: p(Y|A,S,C) = N(Y|A*S,C)
%
% Priors p(A) and p(S) can be arbitrary, because the user must specify
% the posterior estimation function.
%
% Covariance could be a Kronecker product of two arbitrary (or diagonal?)
% covariance matrices C_a and C_s? Missing values cause a nasty problem
% to Kronecker-covariance.. :(
%
% What if also an arbitrary diagonal covariance matrix C_d? The full
% joint covariance matrix would be:
%
% C = sqrt(C_d) * kron(C_a,C_s) * sqrt(C_d)
%
% Then missing values could be considered by using large (infinite)
% values in C_d..?

% Last modified 2010-10-12
% Copyright (c) Jaakko Luttinen (jaakko.luttinen@tkk.fi)

function Q = lcm(Y, D, varargin)

[M,N] = size(Y);

options = struct( ...
    'init',         [],               ...
    'prior',        [],               ...
    'common_tau',   true,             ...
    'common_nu',    false,            ...
    'rotate',       1,                ...
    'update_alpha', 1,               ...
    'update_beta',  1,               ...
    'update_nu',    1,               ...
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

% Helpful variables
Id = speye(D);
log2pi = log(2*pi);


logP = -inf;
loglikelihood = nan*zeros(options.maxiter,1);
start_time = cputime;
time = nan*zeros(options.maxiter,1);

lastsave = now;


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
  [X, XX] = 
  for n=1:N
    Mmv = Obs(:,n);
% $$$     UTau = U(Mmv,n).*tau(Mmv);
% $$$     UWWTau = sum(bsxfun(@times, ...
% $$$                         reshape(U(Mmv,n),[1,1,sum(Mmv)]), ...
% $$$                         WtauW(:,:,Mmv)), ...
% $$$                  3);

    % Distribution
    UWtauW = wsum(WtauW(:,:,Mmv),U(Mmv,n),3);
    CovX(:,:,n) = inv( speye(D) + UWtauW );
    X(:,n) = CovX(:,:,n) * ...
             ( tauW(Mmv,:)'*spdiag(Y(Mmv,n)) - WtauMu(Mmv,:)' ) * U(Mmv,n);

    % Expectations
    XX(:,:,n) = CovX(:,:,n) + X(:,n)*X(:,n)';
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
  % W BLOCK
  %
  %%%%%%%%%%%%%%%%%
  
  %
  % Update distributions and evaluate statistics
  %

  % Update q(W,Mu|tau)
  for m=1:M
    Nmv = Obs(m,:);

    UXX = wsum(XX(:,:,Nmv),U(m,Nmv),3);
    UX = X(:,Nmv)*U(m,Nmv)';
    UY = U(m,Nmv).*Y(m,Nmv);
    UYX = [UY*X(:,Nmv)', sum(UY)];

    % Distribution
    CovWMu(:,:,m) = inv(spdiag([alpha; beta]) + [UXX, UX; UX', sum(U(m,Nmv))]);
    W(m,:) = UYX*CovWMu(:,1:D,m);
    Mu(m) = UYX*CovWMu(:,end,m);
  end
  
  % Update q(tau)
  for m=1:M
    Nmv = Obs(m,:);
    a_tau(m) = p.a_tau + 0.5 * sum(Nmv);
    b_tau(m) = p.b_tau ...
        + 0.5 * (U(m,Nmv).*Y(m,Nmv))*Y(m,Nmv)' ...
        - 0.5 * [W(m,:),Mu(m)] * inv(CovWMu(:,:,m)) * [W(m,:),Mu(m)]';
  end
  if options.common_tau
    % Common tau
    a_tau(:) = p.a_tau + sum(a_tau-p.a_tau);
    b_tau(:) = p.b_tau + sum(b_tau-p.b_tau);
  end
  tau = a_tau./b_tau;
  logtau = psi(a_tau) - log(b_tau);
  for m=1:M
    % Expectations
    tauW(m,:) = tau(m) * W(m,:);
    tauMu(m) = tau(m) * Mu(m);
    WtauW(:,:,m) = tau(m)*W(m,:)'*W(m,:) + CovWMu(1:D,1:D,m);
    MutauMu(m) = tau(m)*Mu(m)^2 + CovWMu(end,end,m);
    WtauMu(m,:) = tau(m)*Mu(m)*W(m,:) + CovWMu(end,1:D,m);
  end
  
  % Update q(alpha)
  if (numel(options.update_alpha)==1 && iter >= options.update_alpha) || ...
        any(iter==options.update_alpha)
    a_alpha(:) = p.a_alpha + 0.5*M;
    b_alpha(:) = p.b_alpha + 0.5*diag(sum(WtauW, 3));
  end
  alpha = a_alpha ./ b_alpha;
  logalpha = psi(a_alpha) - log(b_alpha);
  
  % Update q(beta)
  if (numel(options.update_beta)==1 && iter >= options.update_beta) || ...
        any(iter==options.update_beta)
    a_beta = p.a_beta + 0.5*M;
    b_beta = p.b_beta + 0.5*sum(MutauMu);
  end
  beta = a_beta ./ b_beta;
  logbeta = psi(a_beta) - log(b_beta);
  
  %
  % Compute KL terms
  %
  
  % KL-term: <log q(W,Mu|tau)> - <log p(W|tau,alpha)> - <log p(Mu|tau,beta)>
  KL_WMu = 0;
  for m=1:M
    KL_WMu = KL_WMu - gaussian_entropy(-(D+1)*logtau(m) ...
                                       + logdet_cov(CovWMu(:,:,m)), ...
                                       D);
    KL_WMu = KL_WMu - gaussian_logpdf(diag(WtauW(:,:,m))'*alpha, ...
                                      0, ...
                                      0, ...
                                      -D*logtau(m) - sum(logalpha), ...
                                      D);
    KL_WMu = KL_WMu - gaussian_logpdf(diag(MutauMu(m))'*beta, ...
                                      0, ...
                                      0, ...
                                      -logtau(m) - sum(logbeta), ...
                                      1);
  end
% $$$   KL_W = 0;
% $$$   for m=1:M
% $$$     KL_W = KL_W - gaussian_entropy(-D*logtau(m) + logdet_cov(CovW(:,:,m)), D);
% $$$     KL_W = KL_W - gaussian_logpdf(diag(WtauW(:,:,m))'*w, ...
% $$$                                   0, ...
% $$$                                   0, ...
% $$$                                   -D*logtau(m) - sum(logw), ...
% $$$                                   D);
% $$$   end
  
  % KL-term: <log q(tau)> - <log p(tau)>
  KL_tau = 0;
  if options.common_tau
    KL_tau = KL_tau - gamma_entropy(a_tau(1), b_tau(1));
    KL_tau = KL_tau - gamma_logpdf(tau(1),logtau(1),p.a_tau,p.b_tau);
  else
    KL_tau = KL_tau - sum(gamma_entropy(a_tau, b_tau));
    KL_tau = KL_tau - sum(gamma_logpdf(tau,logtau,p.a_tau,p.b_tau));
  end
  
  % KL-term: <log q(alpha)> - <log p(alpha)>
  KL_alpha = 0;
  KL_alpha = KL_alpha - sum(gamma_entropy(a_alpha, b_alpha));
  KL_alpha = KL_alpha - sum(gamma_logpdf(alpha, logalpha, p.a_alpha, p.b_alpha));
  
  % KL-term: <log q(beta)> - <log p(beta)>
  KL_beta = 0;
  KL_beta = KL_beta - sum(gamma_entropy(a_beta, b_beta));
  KL_beta = KL_beta - sum(gamma_logpdf(beta, logbeta, p.a_beta, p.b_beta));
  
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
      
      E2(m,Nmv) = tau(m)*Y(m,Nmv).^2 ...
          + mtraceprod(WtauW(:,:,m), XX(:,:,Nmv)) ...
          + MutauMu(m) ...
          - 2*Y(m,Nmv).*(tauW(m,:)*X(:,Nmv)) ...
          - 2*Y(m,Nmv).*tauMu(m) ...
          + 2*WtauMu(m,:)*X(:,Nmv);
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
      
      E2(m,Nmv) = tau(m)*Y(m,Nmv).^2 ...
          + mtraceprod(WtauW(:,:,m), XX(:,:,Nmv)) ...
          + MutauMu(m) ...
          - 2*Y(m,Nmv).*(tauW(m,:)*X(:,Nmv)) ...
          - 2*Y(m,Nmv).*tauMu(m) ...
          + 2*WtauMu(m,:)*X(:,Nmv);
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
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% Calculate lower bound of the log-likelihood
  
  monitor_loglikelihood = true;
  
  if monitor_loglikelihood
    
    % Cost from Y
    
    % We can take advantage of the linearity of gaussian_logpdf by summing
    % before the function call:
    tauYh = bsxfun(@plus, tauW*X, tauMu);
    y_invCov_y = (tau(ObsM).*U(Obs))'*(Y(Obs).^2);
    y_invCov_mu = U(Obs)'*(Y(Obs).*tauYh(Obs));
    mu_invCov_mu = 0;
    for m=1:M
      Nmv = Obs(m,:);
      mu_invCov_mu = mu_invCov_mu ...
          + traceprod(WtauW(:,:,m),wsum(XX(:,:,Nmv),U(m,Nmv),3)) ...
          + sum(U(m,Nmv)) * MutauMu(m) ...
          + 2*WtauMu(m,:)*(X(:,Nmv)*U(m,Nmv)');
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
    logP = logpdf_Y - KL_U - KL_X - KL_tau - KL_WMu - KL_alpha - KL_beta;

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
    save(options.autosavefile, 'W','CovW','w','a_w','b_w','X','CovX', ...
         'tau', 'a_tau','b_tau','U','a_U','b_U','nu', 'loglikelihood', ...
         'time', 'iter', 'options');
    lastsave = now;
    fprintf(' done.\n');
  end

  % Transformations for speeding up
  if (numel(options.rotate)==1 && iter>=options.rotate && options.rotate>=1) ...
        || any(iter==options.rotate)
    %orthogonalize();
    rotation2();
  end

end

% $$$ figure
% $$$ kl = [kl_Y(:), kl_U(:), kl_X(:), kl_W(:), kl_tau(:), kl_w(:)];
% $$$ kl = kl(6:end,:);
% $$$ kl = eminus(kl,mean(kl,1));
% $$$ plot(kl)

% Results as a struct
Q.W = W;
Q.CovWMu = CovWMu;
Q.alpha = alpha;
Q.a_alpha = a_alpha;
Q.b_alpha = b_alpha;
Q.Mu = Mu;
Q.beta = beta;
Q.a_beta = beta;
Q.b_beta = b_beta;
Q.X = X;
Q.CovX = CovX;
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

function rotation2()

%alpha_before = alpha

% Rotates W and X, and moves bias from X to Mu.

%Y0 = bsxfun(@plus, W*X, Mu);

sumX = sum(X,2);
sumXX = [sum(XX,3), sumX; sumX', N]; % extended model: [X;1]*[X;1]'
sumWtauW = sum(WtauW,3);
sumMutauMu = sum(MutauMu);
sumWtauMu = sum(WtauMu,1);
sumWMutauWMu = [sumWtauW, sumWtauMu'; sumWtauMu, sumMutauMu];

%diag_wtauw = diag(sum(WtauW,3))

%
% Optimize the rotation
%

R0 = eye([D,D+1]);
%mycheckgrad(@cost, R0(:), 1e-6);
r = minimize(R0(:), @cost, 50);
%r(:) = randn(size(r)); % DEBUG !!!
R = [reshape(r,[D,D+1]); zeros(1,D), 1];
%R
%R = eye(D+1);% DEBUG!!!
invR = inv(R);

%
% Apply the rotation
%

% W and Mu
tauWMu = [tauW,tauMu]*R;
tauW = tauWMu(:,1:D);
tauMu = tauWMu(:,end);
for m=1:M
  Z = [WtauW(:,:,m), WtauMu(m,:)'; WtauMu(m,:), MutauMu(m)];
  % Expectations
  WtauW(:,:,m) = R(:,1:D)'*Z*R(:,1:D);
  MutauMu(m) = R(:,end)'*Z*R(:,end);
  WtauMu(m,:) = R(:,end)'*Z*R(:,1:D);
end


% X
X = invR(1:D,:)*[X;ones(1,N)];
for n=1:N
  XX(:,:,n) = invR(1:D,:)*[XX(:,:,n),X(:,n);X(:,n)',1]*invR(1:D,:)';
end

% alpha and beta
%diag_wtauw = diag(sum(WtauW,3))
b_alpha = p.b_alpha + 0.5*diag(sum(WtauW,3));
alpha = a_alpha ./ b_alpha;
b_beta = p.b_beta + 0.5*sum(MutauMu);
beta = a_beta ./ b_beta;

%alpha_after = alpha

%meanX = mean(X,2)
%sumXX = sum(XX,3) / N
%sumWW = sum(WtauW,3) / M

%Y1 = bsxfun(@plus, W*X, Mu);

%difference = sqrt(mean((Y1(:)-Y0(:)).^2))

%error('jou')


  function [f,df] = cost(r)

  R = [reshape(r,[D,D+1]); zeros(1,D), 1];
  R(end) = 1;
  
  invR = inv(R);
  invRR = inv(R*R');
  logdetR = logdet(R);

  % Transform W and Mu
  WW_R = sumWMutauWMu*R;
  R_WW_R = R'*WW_R;
  
  % Transform q(alpha) and q(beta)  
  % (concatenate alpha and beta to one vector for simplicity)
  z = diag(R_WW_R);
  r_b_alpha = [repmat(p.b_alpha,D,1);p.b_beta] + 0.5*z;
  r_alpha = [a_alpha;a_beta] ./ r_b_alpha;
  r_logalpha = - log(r_b_alpha); % up to a constant
  grad_b = WW_R;
  grad_logb = grad_b * diag(1./r_b_alpha);
  grad_logalpha = -grad_logb;
  grad_alpha = grad_logalpha * diag(r_alpha);
  
  % Transform X

  % Cost from W and Mu
  %
  % TODO: There is something wrong with this!!
  logdet_Cov = 2*M*logdetR;
  x_invCov0_x = traceprod(spdiag(r_alpha),R_WW_R);
  x_invCov0_mu = 0;
  mu_invCov0_mu = 0;
  logdet_Cov0 = -M * sum(r_logalpha);
  grad_logdet_Cov = 2*M*invR';
  grad_x_invCov0_x = 2*WW_R*spdiag(r_alpha) ...
      + grad_alpha * diag(diag(R_WW_R));
  grad_x_invCov0_mu = 0;
  grad_mu_invCov0_mu = 0;
  grad_logdet_Cov0 = -M * grad_logalpha;
  [KL_W,dKL_W] = vb_rotationcost_gaussian(logdet_Cov, ...
                                          x_invCov0_x, ...
                                          x_invCov0_mu, ...
                                          mu_invCov0_mu, ...
                                          logdet_Cov0, ...
                                          grad_logdet_Cov, ...
                                          grad_x_invCov0_x, ...
                                          grad_x_invCov0_mu, ...
                                          grad_mu_invCov0_mu, ...
                                          grad_logdet_Cov0);
  %KL_W = 0;
  %dKL_W = 0;
  
  
  % Cost from alpha and beta
  prior_a = [repmat(p.a_alpha,D,1);p.a_beta];
  prior_b = [repmat(p.b_alpha,D,1);p.b_beta];
  post_a = [a_alpha;a_beta];
  A0_logB = prior_a'*log(r_b_alpha);
  A_B0_invB = prior_b'*r_alpha;
  grad_A0_logB = grad_logb * diag(prior_a);
  grad_A_B0_invB = grad_logalpha * diag(prior_b);
  [KL_alpha,dKL_alpha] = vb_rotationcost_gamma(A0_logB,      ...
                                               A_B0_invB,    ...
                                               grad_A0_logB, ...
                                               grad_A_B0_invB);
  %KL_alpha = 0;
  %dKL_alpha = 0;

  
  % Cost from X
  logdet_Cov = -2*N*logdetR;
  x_invCov0_x = traceprod(invR(1:D,:)'*invR(1:D,:),sumXX);
  x_invCov0_mu = 0;
  mu_invCov0_mu = 0;
  logdet_Cov0 = 0;
  grad_logdet_Cov = -2*N*invR';
  grad_x_invCov0_x = -2*invRR*sumXX*invR';
  grad_x_invCov0_mu = 0;
  grad_mu_invCov0_mu = 0;
  grad_logdet_Cov0 = 0;
  [KL_X,dKL_X] = vb_rotationcost_gaussian(logdet_Cov, ...
                                          x_invCov0_x, ...
                                          x_invCov0_mu, ...
                                          mu_invCov0_mu, ...
                                          logdet_Cov0, ...
                                          grad_logdet_Cov, ...
                                          grad_x_invCov0_x, ...
                                          grad_x_invCov0_mu, ...
                                          grad_mu_invCov0_mu, ...
                                          grad_logdet_Cov0);
  %KL_X = 0;
  %dKL_X = 0;
  
  % Vectorize results
  f = KL_W + KL_alpha + KL_X;
  df = dKL_W + dKL_alpha + dKL_X;
  df = df(1:D,:);
  df = df(:);
  
  
  end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function orthogonalize

  %% Move bias (approximate: mean of X to zero)
% $$$   dmu = mean(X,2);
% $$$   X = X - repmat(dmu,1,N);
% $$$   tauMu = tauMu + tauW*dmu;
% $$$   WtauMu = WtauMu ???
% $$$   MutauMu = ???
    
  %% Rotate W and X
  
  Yh0 = tauW*X;

  % Whitening of X
  [Vx,D2x] = svd(sum(XX,3)/N);
  Dx = spdiag(sqrt(diag(D2x))); % convert to sparse diagonal form
  Qx = Dx \ Vx';
  Qw = full(Vx * Dx); % funny: could be sparse if size(Dx)=[1 1] ??
  tauW = tauW * Qw;
  WtauMu = WtauMu * Qw;
  for m=1:M
%    CovW(:,:,m) = Qw' * CovW(:,:,m) * Qw;
    WtauW(:,:,m) = Qw' * WtauW(:,:,m) * Qw;
% $$$     WW(:,:,m) = W(m,:)'*W(m,:) + CovW(:,:,m);
  end

  % Orthogonalization of W 
  [Vw,Dw] = svd(sum(WtauW,3)/M);
% $$$   [Vw,Dw] = svd(sum(WW,3)/M);
  Qx = Vw' * Qx;
  Qw = Vw;
%  W = W * Qw;
  tauW = tauW * Qw;
  WtauMu = WtauMu * Qw;
  for m=1:M
%    CovW(:,:,m) = Qw' * CovW(:,:,m) * Qw;
    WtauW(:,:,m) = Qw' * WtauW(:,:,m) * Qw; %W(m,:)'*W(m,:) + CovW(:,:,m);
% $$$     WW(:,:,m) = W(m,:)'*W(m,:) + CovW(:,:,m);
  end
  
  % Apply rotations to X
  X = Qx * X;
  for n = 1:N
    CovX(:,:,n) = Qx * CovX(:,:,n) * Qx';
    XX(:,:,n) = X(:,n)*X(:,n)' + CovX(:,:,n);
  end
  
  Yh1 = tauW*X;
  
% $$$   diff_in_orth = norm(Yh1-Yh0)
% $$$   
% $$$   sumWtauW = sum(WtauW,3)
% $$$   sumXX = sum(XX,3) / N
  
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%

function initialize()

  disp('Initialize variables..')

  % Initialize X
  X = nan * zeros(D,N);
  CovX = nan*zeros([D,D,N]); % repmat(eye(D),[1,1,N]);

  % Initialize q(alpha)
  a_alpha = nan * ones(D,1);
  b_alpha = nan * ones(D,1);
  alpha = 1e-3 * ones(D,1);
  logalpha = nan * ones(D,1);
  
  % Initialize q(beta)
  a_beta = nan;
  b_beta = nan;
  beta = 1e-3;
  logbeta = nan;

  % Initialize W and Mu
  W = randn(M,D);
  Mu = zeros(M,1);
  CovWMu = repmat(eye(D+1),[1,1,M]);

  % Initialize tau
  a_tau = nan * ones(M,1);
  b_tau = nan * ones(M,1);
  tau = 1e1 * ones(M,1);
  logtau = nan * ones(M,1);
  
  % Initialize nu and q(U)
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
    % X
    if isfield(options.init, 'X')
      disp('Custom initialization for X.');
      warning(['Currently, X is updated first, so given initialization has ' ...
               'no effect.']);
      X(:,:) = options.init.X;
    end
    if isfield(options.init, 'CovX')
      disp('Custom initialization for CovX.');
      warning(['Currently, X is updated first, so given initialization has ' ...
               'no effect.']);
      CovX(:,:,:) = options.init.CovX;
    end
    % W and Mu
    if isfield(options.init, 'W')
      disp('Custom initialization for W.');
      W(:,:) = options.init.W;
    end
    if isfield(options.init, 'Mu')
      disp('Custom initialization for Mu.');
      Mu(:,:) = options.init.Mu;
    end
    if isfield(options.init, 'CovWMu')
      disp('Custom initialization for CovWMu.');
      CovWMu(:,:,:) = options.init.CovWMu;
    end
    % alpha
    if isfield(options.init, 'alpha')
      disp('Custom initialization for alpha.');
      alpha(:) = options.init.alpha;
    end
% $$$     if isfield(options.init, 'a_alpha')
% $$$       disp('Custom initialization for a_alpha.');
% $$$       a_alpha(:) = options.init.a_alpha;
% $$$     end
% $$$     if isfield(options.init, 'b_alpha')
% $$$       disp('Custom initialization for b_alpha.');
% $$$       b_alpha(:) = options.init.b_alpha;
% $$$     end
    % beta
    if isfield(options.init, 'beta')
      disp('Custom initialization for beta.');
      beta(:) = options.init.beta;
    end
% $$$     if isfield(options.init, 'a_beta')
% $$$       disp('Custom initialization for a_beta.');
% $$$       a_beta(:) = options.init.a_beta;
% $$$     end
% $$$     if isfield(options.init, 'b_beta')
% $$$       disp('Custom initialization for b_beta.');
% $$$       b_beta(:) = options.init.b_beta;
% $$$     end
    % tau
    if isfield(options.init, 'tau')
      disp('Custom initialization for tau.');
      tau(:) = options.init.tau;
    end
% $$$     if isfield(options.init, 'a_tau')
% $$$       disp('Custom initialization for a_tau.');
% $$$       a_tau(:) = options.init.a_tau;
% $$$     end
% $$$     if isfield(options.init, 'b_tau')
% $$$       disp('Custom initialization for b_tau.');
% $$$       b_tau(:) = options.init.b_tau;
% $$$     end
    % nu
    if isfield(options.init, 'nu')
      disp('Custom initialization for nu.');
      nu(:) = options.init.nu;
    end
    % U
    if isfield(options.init, 'U')
      disp('Custom initialization for U.');
      U(:,:) = options.init.U;
    end
  end

  % Initialize second moments for X
  XX = zeros(D,D,N);
  for n=1:N
    XX(:,:,n) = X(:,n)*X(:,n)' + CovX(:,:,n);
  end

  % Initialize expectations for W and Mu
  tauW = zeros(M,D);
  tauMu = zeros(M,1);
  WtauW = zeros(D,D,M);
  MutauMu = zeros(M,1);
  WtauMu = zeros(M,D);
  for m=1:M
    tauW(m,:) = tau(m) * W(m,:);
    tauMu(m) = tau(m) * Mu(m);
    WtauW(:,:,m) = tau(m)*W(m,:)'*W(m,:) + CovWMu(1:D,1:D,m);
    MutauMu(m) = tau(m)*Mu(m)^2 + CovWMu(end,end,m);
    WtauMu(m,:) = tau(m)*Mu(m)*W(m,:) + CovWMu(end,1:D,m);
  end
  
end

%%%%%%%%%%%%%%%%%%%%%%%%%

function set_priors()

  % Default prior
  p.a_alpha = 1e-5;
  p.b_alpha = 1e-5;
  p.a_beta = 1e-5;
  p.b_beta = 1e-5;
  p.a_tau = 1e-5;
  p.b_tau = 1e-5;
  
  % Custom prior
  if isstruct(options.prior)

    % Prior for alpha
    if isfield(options.prior, 'a_alpha')
      p.a_alpha(1) = options.prior.a_alpha;
      disp('Custom prior parameter a_alpha');
    end
    if isfield(options.prior, 'b_alpha')
      p.b_alpha(1) = options.prior.b_alpha;
      disp('Custom prior parameter b_alpha');
    end
    
    % Prior for beta
    if isfield(options.prior, 'a_beta')
      p.a_beta(1) = options.prior.a_beta;
      disp('Custom prior parameter a_beta');
    end
    if isfield(options.prior, 'b_beta')
      p.b_beta(1) = options.prior.b_beta;
      disp('Custom prior parameter b_beta');
    end
    
    % Prior for tau
    if isfield(options.prior, 'a_tau')
      p.a_tau(1) = options.prior.a_tau
      disp('Custom prior parameter a_tau');
    end
    if isfield(options.prior, 'b_tau')
      p.b_tau(1) = options.prior.b_tau
      disp('Custom prior parameter b_tau');
    end
    
  end
  
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function rotate

% I am not absolutely sure that this is bug free

% $$$ D = size(XX,1);
% $$$ N = size(XX,3);
% $$$ M = size(CC,3);

% R should be non-singular because we need to take log(abs(det(R)))

% X_X_ = sum(XX(:,:,1:(N-1)),3) + x0x0;
% CovA = 1/D * (sum(AA,3) - A'*A);
% $$$ Psi = sum(XX,3) - 2*A*sum(X_X,3) + A*X_X_*A' + traceprod(CovA,X_X_)*eye(D);

% Start with identity matrix which does not make the results worse
R = eye(D);

r = minimize(R(:),@cost,100);
R = reshape(r,[D,D]);

% Compute new statistics
for n=1:N
  XX(:,:,n) = R \ XX(:,:,n) / R';
  X_X(:,:,n) = R \ X_X(:,:,n) / R';
end
x0x0 = R \ x0x0 / R';
A = R \ A * R;
CovA = 1/D * trace(inv(R*R')) * R'*CovA*R;
for d=1:D
  AA(:,:,d) = A(d,:)'*A(d,:) + CovA;
end
C = C*R;
for m=1:M
  CC(:,:,m) = R'*CC(:,:,m)*R;
end

I = eye(D);
Rtmp = randn(D);

% Compute new hyperparameters here or after this function!

  function [f,df] = cost(x)
  
    R = reshape(x,[D,D]);
    [f,df] = bound(R);
    f = -f;
    df = -df(:);
    
  end
  
  function [l,dl] = bound(R)

  
  % Use QR to evaluate these more efficiently?
  RR = R'*R;
  invR = inv(R); % inv??
  invRR = inv(R*R'); % inv??
  logdetR = logdet(R);
  trinvRR = trace(invRR);

  Psi = sum(XX,3)*invR'/p.Cov_X - 2*inv(p.Cov_X)*sum(p.mu_X'*X,2);

  % A_invRR_A = A'*invRR*A + trinvRR*CovA;
  % diag_RA_invRR_AR = diag(R'*A_invRR_A*R);
  diagRCCR = diag(R'*sum(WtauW,3)*R); % not efficient

  % Bound terms
  lpX = -0.5 * traceprod(invRR,Psi);
  lqX = N * logdetR;
  lpA = 0; %-0.5 * D * sum(log(diag_RA_invRR_AR));
  lqA = 0; %-0.5*D*logabsdet(trinvRR*RR);
  lpC = -0.5 * M * sum(log(diagRCCR));
  lqC = -M * logdetR;

  % Bound
  l = lpX + lpA + lpC - lqX - lqA - lqC;
  
  
  % Derivative terms
  dlpX = invR * Psi * invRR;
  dlqX = N * invR;
  dlpA = 0; %-D*(diag(1./diag_RA_invRR_AR)*R'*A_invRR_A) ...
         %+ D*invR*A*R*diag(1./diag_RA_invRR_AR)*R'*A'*invRR ...
         %+ D*traceprod(diag(1./diag_RA_invRR_AR),R'*CovA*R)*invR*invRR;
  dlqA = 0; %-D * invR + D^2 * invR*invRR / trinvRR;
  dlpC = -M * (diag(1./diagRCCR)*R'*sum(CC,3));
  dlqC = -M * invR;
  
  % Derivative
  dl = (dlpX + dlpA + dlpC - dlqX - dlqA - dlqC)';
  end

end

end