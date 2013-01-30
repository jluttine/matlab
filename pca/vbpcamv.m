function Q = vbpcamv(Y, d, varargin)
% Modified RPPCA by Jaakko Luttinen
% Uses separate outlier models for each observed signal.
% p(y|W,x,mu,u,tau) = N(y | W*x+mu, diag(1/u*1/tau))
% Uses hierarchical models and variational learning.
%
% If you have troubles with bad minima, you may want to try changing the
% time when the algorithm starts to update the hyperparameters
% (startupdatehyper) and/or to rotate (startrotate).

% Last updated 13th February 2009, Jaakko Luttinen

opts = struct( ...
    'init',          [],...
    'prior', [], ...
    'hiermean', false, ...
    'rotate', true, ...
    'commonnoise', true, ...
    'startrotate', 1, ...
    'startupdatehyper', 1, ...
    'autosavetime', 0,...
    'autosavefile', 'vbpcamv_autosave',...
    'testset', [], ...
    'fixw', false, ...
    'maxiters', 1000, ...
    'convergence', eps);

[ opts, errmsg, wrnmsg ] = argschk( opts, varargin{:} );
if ~isempty(errmsg), error( errmsg ), end
if ~isempty(wrnmsg), warning( wrnmsg ), end

[m,n] = size(Y);

% Missing values
if issparse(Y)
  Obs = Y~=0; %spones(Y);
else
  Obs = ~isnan(Y);
end
[ObsI, ObsJ] = find(Obs);

%%%%%%%%%%%%%%%%%%%
%% INITIALIZATION

nmv = sum(Obs,2);

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
% $$$ mumumu = 0;     % prior for hierarchical stuff
% $$$ taumumu = 1e-5; % prior for hierarchical stuff
v_mumu = 0;
ataumu = 1e-3;  % prior for hierarchical stuff
btaumu = 1e-3;  % prior for hierarchical stuff
mumu = 0;       % init posterior or prior if no hierarchical mean
taumu = 1e-6;   % init posterior or prior if no hierarchical mean
v_mu = 1e0*ones(m,1); % init posterior
mu = zeros(m,1);      % init posterior
  
% Loading matrix (W)
%aw = 1e-16; % prior
%bw = 1e-13; % prior
aw = 1e-10; % prior
bw = 1e-5; % prior
a_w = aw * ones(1,d); % init size
b_w = bw * ones(1,d); % init size
w = 1e-6 * ones(1,d);    % init posterior
logw = nan * ones(1,d);  % init size
W = orth(randn(m,d));          % init posterior
CovW = 1e-0*repmat(eye(d),[1,1,m]); % init posterior
varW = nan*ones(m,d);              % init size

% Observation noise
atau = 1e-4; % prior
btau = 1e-4; % prior
a_tau = nan * ones(m,1);
b_tau = nan * ones(m,1);
tau = 1e2 * ones(m,1); % init posterior (THIS SEEMS TO BE IMPORTANT!!!)
%tau = 1./varmv(Y,2) * (m^2)/((m-d)^2)   % adhoc init posterior
logtau = nan * ones(m,1); % init size
  
% Outlier parameters
% $$$ nu = 10 * ones(m,1); % variable init
% $$$ a_U = nan * ones(m,n); % init size
% $$$ b_U = nan * ones(m,n); % init size
U = ones(m,n);     % init posterior
% $$$ logU = nan * zeros(m,n); % init size

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
    mu(:) = opts.init.mu;
    fprintf('Using old mu.\n');
  end
  if isfield(opts.init, 'v_mu')
    v_mu(:) = opts.init.v_mu;
    fprintf('Using old v_mu.\n');
  end
  if isfield(opts.init, 'mumu')
    mumu(1) = opts.init.mumu;
    fprintf('Using old mumu.\n');
  end
  if isfield(opts.init, 'taumu')
    taumu(1) = opts.init.taumu;
    fprintf('Using old taumu.\n');
  end
  if isfield(opts.init, 'tau')
    tau(:) = opts.init.tau;
    fprintf('Using old tau.\n');
  end
  if isfield(opts.init, 'nu')
    nu(:) = opts.init.nu;
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

Im = eye(m);
Id = eye(d);

lastsave = now;

N = 1:n;
vM = 1:m;
m_mv = sum(Obs,1);
n_mv = sum(Obs,2);
nm_mv = sum(Obs(:)); % number of observations

log2pi = log(2*pi);

oldcosts = -inf;

% Monitor the overall time of the iteration
starttime = cputime;
duration = nan*zeros(opts.maxiters,1);
rmse = nan*zeros(opts.maxiters,1);
rmse_test = nan*zeros(opts.maxiters,1);

cost = nan*zeros(opts.maxiters,1);
for k=1:opts.maxiters
  
  itertime = cputime;
  
  % This is used to monitor the rotation of the subspace after each step
  W_old = W;
  
  % Update X
  for j=1:n
    imv = Obs(:,j);
    UTau = diag(U(imv,j).*tau(imv));
    UWWTau = 0;
    for i=vM(imv)
      UWWTau = UWWTau + U(i,j)*tau(i) * WW(:,:,i);
    end
    CovX(:,:,j) = inv(UWWTau + Id);
    X(:,j) = CovX(:,:,j) * W(imv,:)' * UTau * (Y(imv,j)-mu(imv)); 
    XX(:,:,j) = CovX(:,:,j) + X(:,j)*X(:,j)';
    Sv{j} = CovX(:,:,j); % only for backward compatibility..
  end

  for i=1:m
    jmv = Obs(i,:);

    % Update sum_j(<u><xx>)
    UXX(:,:,i) = 0;
    for j=N(jmv) % sum over observed time instances
      UXX(:,:,i) = UXX(:,:,i) + U(i,j)*XX(:,:,j);%(CovX(:,:,j) + X(:,j)*X(:,j)');
    end

    % Update W
    CovW(:,:,i) = inv(UXX(:,:,i)*tau(i) + diag(w));
    W(i,:) = (U(i,jmv).*(Y(i,jmv)-mu(i))*X(:,jmv)') * CovW(:,:,i) * tau(i);
    WW(:,:,i) = W(i,:)'*W(i,:) + CovW(:,:,i);
      
    % Update mu
    v_mu(i) = 1 / (taumu + sum(U(i,jmv))*tau(i));
    mu(i) = v_mu(i) * (taumu*mumu + U(i,jmv)*(Y(i,jmv)-W(i,:)*X(:,jmv))'*tau(i));

  end
  
  % Rotate
  % THINK/TODO: WHAT WOULD BE A GOOD TIME TO START ROTATING??
  % (rotating too early may lead to bad local minima and not rotating
  % makes converging very slow) should one measure how "converged" the
  % loglikelihood is and then decide whether to start rotate?
  if opts.rotate && k>=opts.startrotate
% $$$     before = trace( sum(CovX,3) * sum(CovW,3) )
    [W,CovW,X,Sv,CovX,mu] = orthogonalize(W,CovW,X,Sv,CovX,mu,opts.fixw,tau,Obs);
% $$$     after = trace( sum(CovX,3) * sum(CovW,3) )
%    [W,CovW,X,Sv,CovX,mu] = orthogonalize_TEST(W,CovW,X,Sv,CovX,mu,tau(1)); false
%    [W,CovW,X,CovX,mu] = RotateToPCA(W,CovW,X,CovX,mu);
  end

  
  for i=1:m
    jmv = Obs(i,:); % observed
    
    WW(:,:,i) = W(i,:)'*W(i,:) + CovW(:,:,i);
    varW(i,:) = diag(CovW(:,:,i));
    
    % Update sum_j(<u><xx>)
    UXX(:,:,i) = 0;
    for j=N(jmv) % sum over observed time instances
      UXX(:,:,i) = UXX(:,:,i) + U(i,j)*(CovX(:,:,j) + X(:,j)*X(:,j)');
    end

    % Update tau
    UWXXW = sum(sum( WW(:,:,i) .* UXX(:,:,i) ));
    a_tau(i) = atau + 0.5*sum(jmv);
    b_tau(i) = btau ...
        + 0.5 * (UWXXW + ...
                 sum(-2*W(i,:)*X(:,jmv)*(U(i,jmv).*(Y(i,jmv)-mu(i)))') + ...
                 sum(U(i,jmv).*(Y(i,jmv)-mu(i)).^2) + ...
                 sum(U(i,jmv)*v_mu(i)));
    
  end
  
  % "AD HOC"(?): Use common variance
  if opts.commonnoise
    a_tau(:) = sum(a_tau-atau) + atau;
    b_tau(:) = sum(b_tau-btau) + btau;
  end
  tau = a_tau./b_tau;
  logtau = psi(a_tau) - log(b_tau);
  % TODO: DEBUGGING: USE ML ESTIMATE
% $$$   logtau = log(tau);

  % Update hyperparameters for mu
  if ~opts.hiermean
    logtaumu = log(taumu);
  else
    a_taumu = ataumu + 0.5*m;
    b_taumu = btaumu + 0.5*sum( (mu-mumu).^2 + v_mu + v_mumu);
    taumu = a_taumu ./ b_taumu; % posterior mean of accuracy hyperparameter
    logtaumu = psi(a_taumu) - log(b_taumu);
  end

  % Update hyperparameter for W (you can try to start updating these
  % hyperparameters only after some iteration if some problems seem to
  % appear.)
  if ~opts.fixw
    if k >= opts.startupdatehyper
      a_w(:) = aw + 0.5*m;
      b_w(:) = bw + 0.5*sum(W.^2 + varW, 1);
    end
    w = a_w ./ b_w; % posterior mean of precision hyperparameter
    logw = psi(a_w) - log(b_w);
  else
    logw = log(w);
  end
  
  % Change of the principal subspace (in radians)
  angle = subspace(W,W_old);

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Calculate lower bound of the log-likelihood

  % Cost from Y
  cost_Y = 0;
  for j=1:n
    imv = Obs(:,j);
    for i=vM(imv)
      mean_muY = W(i,:)*X(:,j)+mu(i);
      var_muY = W(i,:)*CovX(:,:,j)*W(i,:)' + X(:,j)'*CovW(:,:,i)*X(:,j) ...
                + sum(sum( CovX(:,:,j).*CovW(:,:,i) )) + v_mu(i);
      tauY = tau(i);
      logtauY = logtau(i);
      cost_Y = cost_Y + cost_gaussian(Y(i,j),[], mean_muY,var_muY, ...
                                      tauY,logtauY);
    end
  end
    
  % Cost from W
  cost_W = 0;
  for i=1:m
    cost_W = cost_W + cost_gaussian(W(i,:)',CovW(:,:,i), 0,0, w',logw');
  end
    
  % Cost from X
  cost_X = 0;
  for j=1:n
    cost_X = cost_X + cost_gaussian(X(:,j),CovX(:,:,j), 0,0, 1,0);
  end
    
  % Cost from mu
  cost_mu = cost_gaussian(mu,diag(v_mu), mumu,v_mumu, taumu,logtaumu);

  % Cost from tau
  if opts.commonnoise
    cost_tau = cost_gamma(tau(1),logtau(1),a_tau(1),b_tau(1),atau,btau);
  else
    cost_tau = sum(cost_gamma(tau,logtau,a_tau,b_tau,atau,btau));
  end
    
  % Cost from w
  if ~opts.fixw
    cost_w = sum(cost_gamma(w,logw,a_w,b_w,aw,bw));
  else
    cost_w = 0;
  end
  
  if opts.hiermean
    cost_taumu = cost_gamma(taumu,logtaumu,a_taumu,b_taumu,ataumu,btaumu);
  end
    
  oldlogP = logP;
    
  % Lower bound of loglikelihood
  logP = cost_Y + cost_W + cost_X + cost_mu + cost_tau + cost_w;
  %logP = cost_Y;
  costs = [cost_Y; cost_W; cost_X; cost_mu; cost_tau; cost_w];
    
  cost(k) = logP;
    
  % Show progress
  time = cputime - itertime;
  if k > 1
    duration(k) = duration(k-1) + time;
  else
    duration(1) = time;
  end
  fprintf('Step %d: loglike=%e, angle=%e  (%f seconds)\n', ...
          k, logP, angle, time);
  
  % Debugging..
  if logP < oldlogP
%    diff_costs = costs - oldcosts
%    error('Log-likelihood bound not improved! Bug in code or likelihood?');
    warning('Log-likelihood bound decreased! Bug in code or likelihood?');
  end
  oldcosts = costs;
  
  % RMSE
  Yh = W*X + repmat(mu, 1,n);
  err = Y-Yh;
  err = err(~isnan(err));
  rmse(k) = sqrt( mean(err(:).^2) );
  if ~isempty(opts.testset)
    err = opts.testset-Yh;
    err = err(~isnan(err));
    rmse_test(k) = sqrt( mean(err(:).^2) );
  end

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
           'mu','v_mu','tau', 'a_tau','b_tau', 'loglikelihood', 'cost', ...
           'duration', 'rmse', 'rmse_test');
    end
    lastsave = now;
    fprintf(' done.\n');
  end
  
% $$$   if angle < 1e-14
% $$$     fprintf('Stopping iteration: only minor change in subspace\n');
% $$$     break
% $$$   end
  improvement = abs((logP - oldlogP) / logP);
  if improvement < opts.convergence && k > opts.startupdatehyper
    fprintf('Stopping iteration: only minor improvement in loglikelihood\n');
    break
  end

end

fprintf('Saving to %s...', opts.autosavefile);
loglikelihood = logP;
if opts.hiermean
  save(opts.autosavefile, 'W','CovW','w','a_w','b_w','X','CovX','Sv', ...
       'mu','v_mu','mumu','v_mumu','taumu','a_taumu','b_taumu','tau', ...
       'a_tau','b_tau','U','a_U','b_U','nu','loglikelihood');
else
  save(opts.autosavefile, 'W','CovW','w','a_w','b_w','X','CovX','Sv', ...
       'mu','v_mu','tau', 'a_tau','b_tau', 'loglikelihood', 'cost', ...
       'duration', 'rmse', 'rmse_test');
end
fprintf(' done.\n');

% Results as a struct if desired
if true %nargout == 1
  results.W = W;
  results.CovW = CovW;
  results.w = w;
  results.a_w = a_w;
  results.b_w = b_w;
  results.X = X;
  results.CovX = CovX;
  results.Sv = Sv;
  results.mu = mu;
  results.v_mu = v_mu;
  if opts.hiermean
    results.mumu = mumu;
    results.v_mumu = v_mumu;
    results.taumu = taumu;
    results.a_taumu = a_taumu;
    results.b_taumu = b_taumu;
  end
  results.tau = tau;
  results.a_tau = a_tau;
  results.b_tau = b_tau;
% $$$   results.U = U;
% $$$   results.a_U = a_U;
% $$$   results.b_U = b_U;
% $$$   results.nu = nu;
  results.loglikelihood = logP;
  results.cost = cost;
  results.time = duration;
  Q = results;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function c = cost_gaussian(x,Covx,mu,Covmu,tau,logtau)
% c = cost_gaussian(x,Covx,mu,Covmu,tau,logtau)
%
% Calculates the difference between <log p(X)> - <log q(X)>
% where expectation is over q(X).
%
% A lower bound for a Gaussian:
% p(X) = N(MU,1/TAU)
% q(X) = N(x,Covx)
% 
% Rest of the parameters are defined as:
% q(MU) = N(mu,Covmu)
% <TAU> = tau
% <log TAU> = logtau
%
% If Covx==[], then the term <log q(X)> is not calculated.
% This is useful when X is, e.g., observations.


c = 0;

% Cost from q-posterior
if ~isempty(Covx)
  entropy = entropy_gaussian(Covx);
  c = entropy;
else
  Covx = 0;
end

% Cost from prior
err2 = ((x-mu).^2 + diag(Covx) + diag(Covmu));
c = c + 0.5 * sum( -log(2*pi) + logtau - tau .* err2 );


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function c = cost_gamma(x,logx,apost,bpost,aprior,bprior)
% A lower bound for a Gamma distributed X:
% p(X) = G(aprior,bprior)
% q(X) = G(apost,bpost)
% <X> = x
% <log X> = logx

% Cost from prior
c = aprior*log(bprior) - gammaln(aprior) + (aprior-1)*logx - bprior*x;
% Cost from q-posterior
c = c + entropy_gamma(apost,bpost);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function H = entropy_gaussian(Cov)
% Cov is covariance matrix
d = size(Cov,1);

log2pi = log(2*pi);
%H = 0;

%for j=1:size(Cov,3)
L = chol(Cov);
H = d/2*log2pi + 0.5*2*logdettri(L) + d/2;
%end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function H = entropy_gamma(a,b)
% a is shape, b is inverse scale

H = a - log(b) + gammaln(a) - (a-1).*psi(a);
%H = sum(H(:));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ A, Av, S, Sv, mu ] = ...
    RotateToPCA( A, Av, S, Sv, mu);

n1 = size(A,1);
n2 = size(S,2);

% TODO: Take into account the prior for Mu, A???
mS = mean(S,2);
dMu = A*mS;
S = S - repmat(mS,1,n2);
mu = mu + dMu;

covS = S*S';
for j = 1:n2
  covS = covS + Sv(:,:,j);
end

R = 1;

% w.r.t. S
covS = covS / n2;
[VS,D] = eig(covS);
RA = VS*sqrt(D);
A = A*RA;
covA = A'*A;
for i = 1:n1
  Av(:,:,i) = RA'*Av(:,:,i)*RA;
  covA = covA + Av(:,:,i);
end
R = diag(1./sqrt(diag(D)))*VS';

% $$$ % w.r.t. A
% $$$ covA = covA / n1;
% $$$ [VA,DA] = eig(covA);
% $$$ [DA,I] = sort( -diag(DA) );
% $$$ DA = -DA;
% $$$ VA = VA(:,I);
% $$$ A = A*VA;
% $$$ for i = 1:n1
% $$$   Av(:,:,i) = VA'*Av(:,:,i)*VA;
% $$$ end
% $$$ R = VA'*R;

S = R*S;
for j = 1:length(Sv)
    Sv(:,:,j) = R*Sv(:,:,j)*R';
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [W,CovW,X,Sv,CovX,mu] = orthogonalize(W,CovW,X,Sv,CovX,mu,fixw,tau,Obs)
[m,c] = size(W);
n = size(X,2);

WW = W'*W + sum(CovW,3);

% $$$ old_recon = W*X + repmat(mu,1,n);

% Use ML zero mean (mu should be updated after this function!!) ..
if false
  disp('Using complicated but exact mean translation');
  tau = reshape(tau,[1,1,m]);
  nom = 0;
  denom = 0;
  I = eye(c);
  for j=1:n
    imv = Obs(:,j);
    tmp = I + sum(bsxfun(@times,tau(1,1,imv),CovW(:,:,imv)),3);
    nom = nom + tmp * X(:,j);
    denom = denom + tmp;
  end
  dmu = denom \ nom;
else
  dmu = mean(X,2);
end
% .. or my analytic zero mean with proper fix to mu:
% $$$  dmu = inv(n*eye(c)+taumu*WW) * (sum(X,2) + taumu*W'*(mumu-mu));

% Move bias
X = X - repmat(dmu,1,n);
mu = mu + W*dmu;
% $$$ nobiasmove = true

Qx = 1;

% Whiten w.r.t. X
XX = X*X' + sum(CovX,3);
if fixw
  [Vx,Dx,tmp] = svd(XX/(n-m)); % USE THIS IF FIXED w ??
else
  [Vx,Dx,tmp] = svd(XX/n);
end
Qx = diag(1./sqrt(diag(Dx))) * Vx';
Qw = Vx*sqrt(Dx);
W = W * Qw;
for i=1:size(CovW,3)
  CovW(:,:,i) = Qw'*CovW(:,:,i)*Qw;
end
X = Qx * X;
for j = 1:size(CovX,3)
  Sv{j} = Qx*Sv{j}*Qx';
  CovX(:,:,j) = Qx*CovX(:,:,j)*Qx';
end

% Check that XX is really whitened! (because of numerical issues!!)
XX = X*X' + sum(CovX,3);
if fixw
  [Vx,Dx,tmp] = svd(XX/(n-m)); % USE THIS IF FIXED w ??
else
  [Vx,Dx,tmp] = svd(XX/n);
end
if Dx(1) > 1.1
  warning('Needs to whiten X again. See vbpcamv->orthogonalize');
  % Whiten w.r.t. X AGAIN
  Qx = diag(1./sqrt(diag(Dx))) * Vx';
  Qw = Vx*sqrt(Dx);
  W = W * Qw;
  for i=1:size(CovW,3)
    CovW(:,:,i) = Qw'*CovW(:,:,i)*Qw;
  end
  X = Qx * X;
  for j = 1:size(CovX,3)
    Sv{j} = Qx*Sv{j}*Qx';
    CovX(:,:,j) = Qx*CovX(:,:,j)*Qx';
  end
end

Qx = 1;
% Diagonalize w.r.t. W
WW = W'*W + sum(CovW,3);
[Vw,Dw,tmp] = svd(WW);
%[Dw,I] = sort(diag(Dw), 'descend');
%Vw = Vw(:,I);
Qx = Vw' * Qx;
Qw = Vw;
W = W * Qw;
for i=1:size(CovW,3)
  CovW(:,:,i) = Qw'*CovW(:,:,i)*Qw;
end
X = Qx * X;
for j = 1:size(CovX,3)
  Sv{j} = Qx*Sv{j}*Qx';
  CovX(:,:,j) = Qx*CovX(:,:,j)*Qx';
end

% $$$ xx = X*X' + sum(CovX,3)
% $$$ ww = W'*W + sum(CovW,3)

% $$$ new_recon = W*X + repmat(mu,1,n);
% $$$ 
% $$$ mse = mean((old_recon(:) - new_recon(:)).^2)

return
% $$$ %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% $$$ function [W,CovW,X,Sv,CovX,mu] = orthogonalize_TEST(W,CovW,X,Sv,CovX,mu,tau)
% $$$ [m,c] = size(W);
% $$$ n = size(X,2);
% $$$ 
% $$$ WW = W'*W + sum(CovW,3);
% $$$ 
% $$$ old_recon = W*X + repmat(mu,1,n);
% $$$ 
% $$$ % Use ML zero mean (mu should be updated after this function!!) ..
% $$$ % $$$ dmu = mean(X,2);
% $$$ % $$$ tau*sum(CovW,3)
% $$$ dmu = inv(n*eye(c)+tau*n*sum(CovW,3)) * (sum((eye(c)+tau*sum(CovW,3))*X, 2));
% $$$ %(dmu2-dmu)./dmu
% $$$ % .. or my analytic zero mean with proper fix to mu:
% $$$ % $$$  dmu = inv(n*eye(c)+taumu*WW) * (sum(X,2) + taumu*W'*(mumu-mu));
% $$$ 
% $$$ % Move bias
% $$$ X = X - repmat(dmu,1,n);
% $$$ mu = mu + W*dmu;
% $$$ 
% $$$ Qx = 1;
% $$$ 
% $$$ % Whiten w.r.t. X
% $$$ XX = X*X' + sum(CovX,3);
% $$$ [Vx,Dx] = eig(XX/n);
% $$$ %[Vx,Dx] = eig(XX/(n+m));
% $$$ %[Vx,Dx] = eig(XX/(n-m));
% $$$ Qx = diag(1./sqrt(diag(Dx))) * Vx';
% $$$ Qw = Vx*sqrt(Dx);
% $$$ W = W * Qw;
% $$$ for i=1:size(CovW,3)
% $$$   CovW(:,:,i) = Qw'*CovW(:,:,i)*Qw;
% $$$ end
% $$$ 
% $$$ % Whiten w.r.t. W
% $$$ WW = W'*W + sum(CovW,3);
% $$$ [Vw,Dw] = eig(WW);
% $$$ [Dw,I] = sort(diag(Dw), 'descend');
% $$$ Vw = Vw(:,I);
% $$$ Qx = Vw' * Qx;
% $$$ Qw = Vw;
% $$$ W = W * Qw;
% $$$ for i=1:size(CovW,3)
% $$$   CovW(:,:,i) = Qw'*CovW(:,:,i)*Qw;
% $$$ end
% $$$ 
% $$$ X = Qx * X;
% $$$ for j = 1:size(CovX,3)
% $$$   Sv{j} = Qx*Sv{j}*Qx';
% $$$   CovX(:,:,j) = Qx*CovX(:,:,j)*Qx';
% $$$ end
% $$$ 
% $$$ new_recon = W*X + repmat(mu,1,n);
% $$$ 
% $$$ mse = mean((old_recon(:) - new_recon(:)).^2);
% $$$ 
% $$$ return
% $$$ 
