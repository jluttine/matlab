function Q = rvbpcamv(Y, d, varargin)
% Modified RPPCA by Jaakko Luttinen
% Uses separate outlier models for each observed signal.
% p(y|W,x,mu,u,tau) = N(y | W*x+mu, diag(1/u*1/tau))
% Uses hierarchical models and variational learning.
%
% If you have troubles with bad minima, you may want to try changing the
% time when the algorithm starts to update the hyperparameters
% (startupdatehyper) and/or to rotate (startrotate).

warning("This method is deprecated. See fa/demo_rfa.m")

opts = struct( ...
    'init',          [],...
    'prior', [], ...
    'hiermean', false, ...
    'rotate', true, ...
    'commonnoise', true, ...
    'common_nu', false, ...
    'startrotate', 1, ...
    'startupdatehyper', 5, ...
...%    'bias',          'update',... % { 'none'|'rowmean' }
...%    'uniquesv',      1,......%    'autosave',      3600,...
...%    'minangle',      1e-8,...
...%    'verbose',       1,...
...%    'display',       0 );
    'autosavetime',      0,...
    'autosavefile', 'rvbpcamv_autosave',...
    'maxiters',      100);

[ opts, errmsg, wrnmsg ] = argschk( opts, varargin{:} );
if ~isempty(errmsg), error( errmsg ), end
if ~isempty(wrnmsg), warning( wrnmsg ), end

[m,n] = size(Y);

% Missing values
M = isnan(Y);
Obs = ~M;
[ObsI, ObsJ] = find(Obs);

%%%%%%%%%%%%%%%%%%%
%% INITIALIZATION

nmv = sum(~M,2);

% NAMING CONVENTION:
% variables are named without '_': e.g. a->tau->mu ==> ataumu
% but the values of posterior distribution parameters are named with one
% '_': e.g. a->tau->mu ==> a_taumu
% So, priors (or hyperparameter variables) have no '_'. Understand the
% difference!
% But for clarity mu_W is just W and mu_X is just X and so on. That is,
% plain variable name (W,X,mu,tau,logtau,U,logU) refers to posterior
% expectation if there is no naming conflict.


% Mean (mu)
mumumu = 0;     % prior for hierarchical stuff
taumumu = 1e-5; % prior for hierarchical stuff
ataumu = 1e-5;  % prior for hierarchical stuff
btaumu = 1e-5;  % prior for hierarchical stuff
mumu = 0;       % init posterior or prior if no hierarchical mean
taumu = 1e-5;   % init posterior or prior if no hierarchical mean
v_mu = 1e5*ones(m,1); % init posterior
mu = zeros(m,1);      % init posterior
  
% Loading matrix (W)
%aw = 1e-16; % prior
%bw = 1e-13; % prior
aw = 1e-6; % prior
bw = 1e-3; % prior
a_w = aw * ones(1,d); % init size
b_w = bw * ones(1,d); % init size
w = 1e-3 * ones(1,d);    % init posterior
logw = nan * ones(1,d);  % init size
W = orth(randn(m,d));          % init posterior
CovW = 1e-0*repmat(eye(d),[1,1,m]); % init posterior
varW = nan*ones(m,d);              % init size

% Observation noise
atau = 1e-5; % prior
btau = 1e-5; % prior
a_tau = nan * ones(m,1);
b_tau = nan * ones(m,1);
tau = 1e2 * ones(m,1); % init posterior (THIS SEEMS TO BE IMPORTANT!!!)
%tau = 1./varmv(Y,2) * (m^2)/((m-d)^2)   % adhoc init posterior
logtau = nan * ones(m,1); % init size
  
% Outlier parameters
nu = 10 * ones(m,1); % variable init
a_U = nan * ones(m,n); % init size
b_U = nan * ones(m,n); % init size
U = ones(m,n);     % init posterior
logU = nan * zeros(m,n); % init size

% Initialize sizes
X = nan * zeros(d,n); % init size
Sv = cell(n,1);
Sv(:) = {nan*zeros(d,d)};
CovX = nan*zeros(d,d,n);
XX = nan*zeros(d,d,n);
UXX = nan*zeros(d,d,m);


% Prior information
if isstruct(opts.prior)
  if isfield(opts.prior, 'mumumu')
    mumumu(1) = opts.prior.mumumu
  end
  if isfield(opts.prior, 'taumumu')
    taumumu(1) = opts.prior.taumumu
  end
  if isfield(opts.prior, 'ataumu')
    ataumu(1) = opts.prior.ataumu
  end
  if isfield(opts.prior, 'btaumu')
    btaumu(1) = opts.prior.btaumu
  end
  if isfield(opts.prior, 'aw')
    aw(1) = opts.prior.aw
  end
  if isfield(opts.prior, 'bw')
    bw(1) = opts.prior.bw
  end
  if isfield(opts.prior, 'atau')
    atau(1) = opts.prior.atau
  end
  if isfield(opts.prior, 'btau')
    btau(1) = opts.prior.btau
  end
end

% Use given initialization (size matching is checked by using ':' )
if isstruct(opts.init)
  if isfield(opts.init, 'W')
    W(:,:) = opts.init.W;
    fprintf('Using old W.\n');
  end
  if isfield(opts.init, 'CovW')
    CovW(:,:,:) = opts.init.CovW;
    fprintf('Using old CovW.\n');
  end
  if isfield(opts.init, 'w')
    w(:) = opts.init.w
    fprintf('Using old w.\n');
  end
  if isfield(opts.init, 'X')
    X(:,:) = opts.init.X;
    fprintf('Using old X.\n');
  end
  if isfield(opts.init, 'CovX')
    CovX(:,:,:) = opts.init.CovX;
    fprintf('Using old CovX.\n');
  end
  if isfield(opts.init, 'Sv')
    Sv(:) = opts.init.Sv;
    fprintf('Using old Sv.\n');
  end
  if isfield(opts.init, 'mu')
    mu(:) = opts.init.mu
    fprintf('Using old mu.\n');
  end
  if isfield(opts.init, 'v_mu')
    v_mu(:) = opts.init.v_mu
    fprintf('Using old v_mu.\n');
  end
  if isfield(opts.init, 'mumu')
    mumu(1) = opts.init.mumu
    fprintf('Using old mumu.\n');
  end
  if isfield(opts.init, 'taumu')
    taumu(1) = opts.init.taumu
    fprintf('Using old taumu.\n');
  end
  if isfield(opts.init, 'tau')
    tau(:) = opts.init.tau
    fprintf('Using old tau.\n');
  end
  if isfield(opts.init, 'nu')
    nu(:) = opts.init.nu
    fprintf('Using old nu.\n');
  end
  if isfield(opts.init, 'U')
    U(:,:) = opts.init.U;
    fprintf('Using old U.\n');
  end
end

WW = zeros(d,d,m);
for i=1:m
  WW(:,:,i) = W(i,:)'*W(i,:) + CovW(:,:,i);
end

% Log-likelihood
logP = -Inf;

Im = speye(m);
Id = speye(d);

lastsave = now;

N = 1:n;
vM = 1:m;
m_mv = sum(Obs,1);
n_mv = sum(Obs,2);
nm_mv = sum(Obs(:)); % number of observations

log2pi = log(2*pi);

% Mean squared errors < (y-wx-m)^2 >:
E2 = zeros(size(Y));

cost = nan*zeros(opts.maxiters,1);
for k=1:opts.maxiters
  
  tic;
  
  % This is used to monitor the rotation of the subspace after each step
  W_old = W;
  
  % Update X
  for j=1:n
    imv = ~M(:,j);
    UTau = U(imv,j).*tau(imv);
    UWWTau = 0;
    for i=vM(imv)
      UWWTau = UWWTau + U(i,j)*tau(i) * WW(:,:,i);
    end
    CovX(:,:,j) = inv(UWWTau + Id);
    X(:,j) = CovX(:,:,j) * ( W(imv,:)' * ( UTau .* (Y(imv,j)-mu(imv)) ) ); 
    XX(:,:,j) = CovX(:,:,j) + X(:,j)*X(:,j)';
    %Sv{j} = CovX(:,:,j); % only for backward compatibility..
  end

  % Update W
  for i=1:m
    jmv = ~M(i,:);

    % Update sum_j(<u><xx>)
    UXX(:,:,i) = 0;
    for j=N(jmv) % sum over observed time instances
      UXX(:,:,i) = UXX(:,:,i) + U(i,j)*XX(:,:,j);%(CovX(:,:,j) + X(:,j)*X(:,j)');
    end

    % Update W
    CovW(:,:,i) = inv(UXX(:,:,i)*tau(i) + diag(w));
    W(i,:) = (U(i,jmv).*(Y(i,jmv)-mu(i))*X(:,jmv)') * CovW(:,:,i) * tau(i);
    WW(:,:,i) = W(i,:)'*W(i,:) + CovW(:,:,i);
      
  end
  
  % Update mu
  for i=1:m
    jmv = ~M(i,:); % observed
    
    WW(:,:,i) = W(i,:)'*W(i,:) + CovW(:,:,i);
    varW(i,:) = diag(CovW(:,:,i));
    
    v_mu(i) = 1 / (taumu + sum(U(i,jmv))*tau(i));
    mu(i) = v_mu(i) * (taumu*mumu + U(i,jmv)*(Y(i,jmv)-W(i,:)*X(:,jmv))'* ...
                       tau(i));
  end

  % Evaluate squared errors
  for i=1:m
    jmv = Obs(i,:);
    for j=N(jmv)
      E2(i,j) = (Y(i,j)-mu(i))^2 + v_mu(i) - 2*(Y(i,j)-mu(i))*(W(i,:)*X(:,j)) ...
                + traceprod(WW(:,:,i),XX(:,:,j));
    end
  end
  
  % Update nu
  if k >= opts.startupdatehyper
    if opts.common_nu
      error('Common nu not implemented yet.');
    else
      % Common degrees of freedom
      for i=1:m
        if false
          disp('Old style nu update')
          % old style
          c = mean(logU(i,jmv)-U(i,jmv)); % separate degrees ignoring m.v.
          nu(i) = exp(fminsearch(@(nu) ((1+log(exp(nu)/2)-psi(exp(nu)/2)+c)^2), ...
                                 log(nu(i))));
        else
          % new style (type-2 ML)
          nu(i) = t_ml(tau(i)*E2(i,~M(i,:)), nu(i));
        end
      end
    end
  end
  

  % Update U
  for i=1:m
    jmv = Obs(i,:);
    a_U(i,jmv) = (nu(i)+1) / 2;
    b_U(i,jmv) = (nu(i) + E2(i,jmv)*tau(i)) / 2;
  end
% $$$   a_U(Obs) = (nu(ObsI)+1) / 2;
% $$$   b_U(Obs) = (nu(ObsI) + E2(Obs).*tau(ObsI)) / 2;
  U(Obs) = a_U(Obs)./b_U(Obs);
  logU(Obs) = psi(a_U(Obs)) - log(b_U(Obs));
  
  % Update tau
  if opts.commonnoise
    a_tau(:) = atau + 0.5 * sum(~M(:));
    b_tau(:) = btau + 0.5 * U(~M(:))'*E2(~M(:));
  else
      
    for i=1:m
      jmv = ~M(i,:);
      
      a_tau(i) = atau + 0.5 * sum(jmv);
      b_tau(i) = btau + 0.5 * U(i,jmv)*E2(i,jmv)';
    end
  end
  tau = a_tau./b_tau;
  logtau = psi(a_tau) - log(b_tau);

  % Rotate
  if opts.rotate && k > opts.startrotate
    orthogonalize;
% $$$     [W,CovW,X,Sv,CovX,mu] = orthogonalize(W,CovW,X,Sv,CovX,mu);
  end

  % Update hyperparameter for W (you can try to start updating these
  % hyperparameters only after some iteration if some problems seem to
  % appear.)
  if k >= opts.startupdatehyper
    a_w(:) = aw + 0.5*m;
    b_w(:) = bw + 0.5*sum(W.^2 + varW, 1);
  end
  w = a_w ./ b_w; % posterior mean of accuracy hyperparameter
  logw = psi(a_w) - log(b_w);
  
  % Change in the update of the principal subspace
  angle = subspace(W,W_old);

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Calculate lower bound of the log-likelihood
  if true
% $$$     err = 0;
% $$$     for i=1:m
% $$$       jmv = ~M(i,:);
% $$$ 
% $$$       UXX(:,:,i) = 0;
% $$$       for j=N(jmv) % sum over observed time instances
% $$$         UXX(:,:,i) = UXX(:,:,i) + U(i,j)*(CovX(:,:,j) + X(:,j)*X(:,j)');
% $$$       end
% $$$ 
% $$$       % Update WW and varW
% $$$       WW(:,:,i) = W(i,:)'*W(i,:) + CovW(:,:,i);
% $$$       varW(i,:) = diag(CovW(:,:,i));
% $$$       
% $$$       UWXXW = sum(sum( WW(:,:,i) .* UXX(:,:,i) ));
% $$$       err = err + tau(i) * ...
% $$$             (UWXXW ...
% $$$              + sum(-2*W(i,:)*X(:,jmv)*(U(i,jmv).*(Y(i,jmv)-mu(i)))') ...
% $$$              + sum(U(i,jmv).*(Y(i,jmv)-mu(i)).^2) ...
% $$$              + sum(U(i,jmv)*v_mu(i)));
% $$$     
% $$$     end
% $$$ 
    
    % Cost from Y
    cost_Y = -0.5*(tau(ObsI).*U(Obs))'*E2(Obs) - 0.5*sum(log2pi - logtau(ObsI) ...
                                                      - logU(Obs));

    % Cost from W
    cost_W = 0;
    for i=1:m
      cost_W = cost_W + mnorm_vbcost(W(i,:)',chol(CovW(:,:,i))', 0,0, ...
                                     spdiag(1./sqrt(w)), -sum(logw));
    end

    % Cost from X
    cost_X = 0;
    for j=1:n
      cost_X = cost_X + mnorm_vbcost(X(:,j),chol(CovX(:,:,j))', 0,0,Id,0);
    end

    % Cost from mu
    cost_mu = 0;
    for i=1:m
      cost_mu = cost_mu + mnorm_vbcost(mu(i),sqrt(v_mu(i)), mumu,0, ...
                                       1./sqrt(taumu), -log(taumu));
    end

    % Cost from tau
    if opts.commonnoise
      cost_tau = gam_vbcost(tau(1),logtau(1),a_tau(1), b_tau(1),atau,btau);
    else
      cost_tau = sum( gam_vbcost(tau,logtau,a_tau, b_tau, atau,btau) );
    end

    % Cost from w
    cost_w = sum( gam_vbcost(w,logw,a_w,b_w,aw,bw) );
    
    % Cost from U
    cost_U = sum( gam_vbcost(U(Obs),logU(Obs), a_U(Obs), b_U(Obs), nu(ObsI)/2, ...
                             nu(ObsI)/2 ) );
    
    
    oldlogP = logP;
    
    logP = cost_Y + cost_W + cost_X + cost_mu + cost_tau + cost_w + cost_U;

    cost(k) = logP;

    % Debugging..
    if logP < oldlogP
      warning('Log-likelihood bound not improved! Bug in code or likelihood?');
    end
    
  end
  
  % Show progress
  time = toc;
  fprintf('Step %d: loglike=%e, angle=%e  (%f seconds)\n', ...
          k, logP, angle, time);
  
  % Check whether to save the results
  if (now-lastsave)*3600*24 >= opts.autosavetime && opts.autosavetime > 0
    fprintf('Saving to %s...', opts.autosavefile);
    loglikelihood = logP;
    if opts.hiermean
      save(opts.autosavefile, 'W','CovW','w','a_w','b_w','X','CovX','Sv', ...
           'mu','v_mu','mumu','v_mumu','taumu','a_taumu','b_taumu','tau', ...
           'a_tau','b_tau','U','a_U','b_U','nu','loglikelihood');
    else
      save(opts.autosavefile, 'W','CovW','w','a_w','b_w','X','CovX','Sv', ...
           'mu','v_mu','tau', 'a_tau','b_tau','U','a_U','b_U','nu', ...
           'loglikelihood');
    end
    lastsave = now;
    fprintf(' done.\n');
  end

% $$$   if angle < 1e-15
% $$$     fprintf('Stopping iteration: only minor change in subspace\n');
% $$$     break
% $$$   end
% $$$   improvement = abs((logP - oldlogP) / logP);
% $$$   if improvement < 1e-8
% $$$     fprintf('Stopping iteration: only minor improvement in loglikelihood\n');
% $$$     break
% $$$   end

end

% Results as a struct if desired
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
Q.loglikelihood = logP;
Q.cost = cost;


%% NESTED FUNCTIONS

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function orthogonalize
%function [W,CovW,X,Sv,CovX,mu] = orthogonalize(W,CovW,X,Sv,CovX,mu)

  %% Move bias (approximate: mean of X to zero)
  dmu = mean(X,2);
  X = X - repmat(dmu,1,n);
  mu = mu + W*dmu;

  %% Rotate W and X

  % Whitening of X
  XX = X*X' + sum(CovX,3);
  [Vx,D2x] = svd(sum(XX,3)/n);
  Dx = spdiag(sqrt(diag(D2x))); % convert to sparse diagonal form
  Qx = Dx \ Vx';
  Qw = full(Vx * Dx); % funny: could be sparse if size(Dx)=[1 1] ??
  W = W * Qw;
  for i=1:m
    CovW(:,:,i) = Qw' * CovW(:,:,i) * Qw;
    WW(:,:,i) = W(i,:)'*W(i,:) + CovW(:,:,i);
  end

  % Orthogonalization of W 
  [Vw,Dw] = svd(sum(WW,3)/m);
  Qx = Vw' * Qx;
  Qw = Vw;
  W = W * Qw;
  for i=1:m
    CovW(:,:,i) = Qw' * CovW(:,:,i) * Qw;
    WW(:,:,i) = W(i,:)'*W(i,:) + CovW(:,:,i);
  end
  
  % Apply rotations to X
  X = Qx * X;
  for j = 1:n
    % Sv{j} = Qx*Sv{j}*Qx';
    CovX(:,:,j) = Qx * CovX(:,:,j) * Qx';
    XX(:,:,j) = X(:,j)*X(:,j)' + CovX(:,:,j);
  end
  
end

% $$$ function orthogonalize
% $$$ %function [W,CovW,X,Sv,CovX,mu] = orthogonalize(W,CovW,X,Sv,CovX,mu)
% $$$ 
% $$$   [m,c] = size(W);
% $$$   n = size(X,2);
% $$$ 
% $$$   %% Move bias (approximate: mean of X to zero)
% $$$   dmu = mean(X,2);
% $$$   X = X - repmat(dmu,1,n);
% $$$   mu = mu + W*dmu;
% $$$ 
% $$$   %% Rotate W and X
% $$$ 
% $$$   %WW = W'*W + sum(CovW,3);
% $$$ 
% $$$   Qx = 1;
% $$$ 
% $$$   % Whiten X
% $$$   XX = X*X' + sum(CovX,3);
% $$$   [Vx,Dx] = eig(XX/n);
% $$$   Qx = inv(sqrt(Dx)) * Vx';
% $$$   Qw = Vx*sqrt(Dx);
% $$$   W = W * Qw;
% $$$   for i=1:size(CovW,3)
% $$$     CovW(:,:,i) = Qw'*CovW(:,:,i)*Qw;
% $$$   end
% $$$ 
% $$$   % Whiten w.r.t. W
% $$$   WW = W'*W + sum(CovW,3);
% $$$   [Vw,Dw] = eig(WW/m);
% $$$   [Dw,I] = sort(diag(Dw), 'descend');
% $$$   Vw = Vw(:,I);
% $$$   Qx = Vw' * Qx;
% $$$   Qw = Vw;
% $$$   W = W * Qw;
% $$$   for i=1:size(CovW,3)
% $$$     CovW(:,:,i) = Qw'*CovW(:,:,i)*Qw;
% $$$   end
% $$$   X = Qx * X;
% $$$   for j = 1:length(Sv)
% $$$     Sv{j} = Qx*Sv{j}*Qx';
% $$$     CovX(:,:,j) = Qx*CovX(:,:,j)*Qx';
% $$$   end
% $$$ 
% $$$ end
% $$$ 

% $$$ %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% $$$ 
% $$$ function H = entropy_gaussian(Cov)
% $$$ % Cov is covariance matrix
% $$$ d = size(Cov,1);
% $$$ 
% $$$ log2pi = log(2*pi);
% $$$ H = 0;
% $$$ 
% $$$ for j=1:size(Cov,3)
% $$$   H = H + d/2*log2pi + 0.5*log(det(Cov(:,:,j))) + d/2;
% $$$ end
% $$$ 
% $$$ return
% $$$ 
% $$$ %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% $$$ 
% $$$ function H = entropy_gamma(a,b)
% $$$ % a is shape, b is inverse scale
% $$$ 
% $$$ H = a - log(b) + gammaln(a) - (a-1).*psi(a);
% $$$ H = sum(H(:));
% $$$ 
% $$$ return

end