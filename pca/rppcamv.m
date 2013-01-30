function [W,X,invB,mu,v,tau,u,a_u,b_u,loglikelihood] = rppcamv(Y, d, varargin)
% Implementation of robust probabilistic PCA with missing values.
% See Archambeau, Delannay, Verleysen: "Robust Probabilistic
% Projections", 2006.
%
% [W,X,Sv,mu,v,tau,u] = rppcamv(Y, d)
%
% Y is the m x n data matrix with missing values marked as NaNs.
% d is the number of dimensions of the latent space.

opts = struct( ...
    'init',          'random',...
    'maxiters',      100,...
    'rotate', true, ...
    'testset', [], ...
...%    'bias',          'update',... % { 'none'|'rowmean' }
...%    'uniquesv',      1,...
    'autosavetime',      0,...
    'autosavefile', 'rppcamv_autosave',...
...%    'minangle',      1e-8,...
...%    'verbose',       1,...
...%    'display',       0 );
    'outliers', true, ...
    'noise', true);

[ opts, errmsg, wrnmsg ] = argschk( opts, varargin{:} );
if ~isempty(errmsg), error( errmsg ), end
if ~isempty(wrnmsg), warning( wrnmsg ), end

[D,n] = size(Y);

% Boolean matrix of missing values
M = isnan(Y);

% initializing guesses
mu = zeros(D,1);%mean(Y,2);
Mu = repmat(mu,1,n);
W = orth(randn(D,d));
tau = 1;
if opts.outliers
  v = 10;
else
  v = 1e100;
end

% initialize sizes
u = zeros(1,n);
logu = zeros(1,n);
a_u = nan * u;
b_u = nan * u;
X = zeros(d,n);
S = zeros(d,d,n);
ID = eye(D);
Id = eye(d);

invB = zeros(d,d,n);

lastsave = now;

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
    w(:) = opts.init.w;
    fprintf('Using old w.\n');
  end
  if isfield(opts.init, 'mu')
    mu(:) = opts.init.mu
    fprintf('Using old mu.\n');
  end
  if isfield(opts.init, 'tau')
    tau(:) = opts.init.tau;
    fprintf('Using old tau.\n');
  end
  if isfield(opts.init, 'v')
    v(:) = opts.init.v;
    fprintf('Using old v.\n');
  end
  if isfield(opts.init, 'u')
    u(:,:) = opts.init.u;
    fprintf('Using old u.\n');
  end
end

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
  
  %%%%%%%%%%
  %% E step
  for j=1:n
    imv = M(:,j);   % zero-one indeces of missing dimensions
    omv = ~imv;     % zero-one indeces of observed dimensions
    Dmv = sum(omv); % number of observed dimensions
    % Update u
    A = inv(W(omv,:)*W(omv,:)' + (1/tau)*eye(Dmv));
    a_u(j) = (Dmv+v)/2;
    b_u(j) = ((Y(omv,j)-mu(omv))'*A*(Y(omv,j)-mu(omv)) + v) / 2;
    u(j) = a_u(j)/b_u(j);
    logu(j) = psi(a_u(j)) - log(b_u(j));
    % Update x
    invB(:,:,j) = inv(tau * W(omv,:)' * W(omv,:) + Id);
    X(:,j) = tau * invB(:,:,j) * W(omv,:)' * (Y(omv,j)-mu(omv)); % E{x}
    S(:,:,j) = invB(:,:,j) + u(j)*X(:,j)*X(:,j)';         % E{uxx'}
  end
  
  %%%%%%%%%%
  %% M STEP 
  %% Update mu
  for i=1:D
    jmv = ~M(i,:);
    if any(jmv)
      % If any observations
      mu(i) = sum(u(jmv).*(Y(i,jmv)-W(i,:)*X(:,jmv))) / sum(u(jmv));
    else
      % If no observations. Hierarchical model could give something to this..
      mu(i) = 0; % ???
    end
  end
  Mu = repmat(mu,1,n);
  %% Update W
  for i=1:D
    jmv = ~M(i,:);
    if any(jmv)
      W(i,:) = u(jmv).*(Y(i,jmv)-mu(i))*X(:,jmv)' * inv(sum(S(:,:,jmv),3));
    else
      % No observations. What to do? Hierarchical model has no problems
      % with this..
      W(i,:) = 0;
    end
  end
  
  %% Update noise
  Nmv = sum(~M(:));
  trSWW = 0;
  midterm = 0;
  for j=1:n
    imv = ~M(:,j);
    trSWW = trSWW + sum(sum( S(:,:,j).*(W(imv,:)'*W(imv,:)) ));
    midterm = midterm - 2*u(j) * (Y(imv,j)-mu(imv))'*W(imv,:)*X(:,j);
  end
  Err = Y - Mu;
  Err(M) = 0;
  err = sum(sum( (Err.^2) * u' ));
  tau = Nmv / (err + midterm + trSWW);

  %% Update v (degrees of freedom)
  if opts.outliers
    c = mean(logu-u);
    v = exp(fminsearch(@(v) ((1+log(exp(v)/2)-psi(exp(v)/2)+c)^2), log(v)));
  end
  
  
  %%%%%%%%%%%%%%%
  %% OTHER STUFF
  
  % Orthogonalize the columns of W
  if opts.rotate
    [X,invB,S,W,mu] = orthogonalize(X,u,invB,S,W,mu);
  end

  % The change (in radians) in the subspace at one step
  angle = subspace(W, W_old);

  % Calculate log-likelihood
  logP = rppcamv_loglikelihood(Y, W, [], [], mu, tau, 0.5*v.*ones(n,1), ...
                               0.5*v.*ones(n,1));
  
  cost(k) = logP;
    
  % Show progress
  time = cputime - itertime;
  if k > 1
    duration(k) = duration(k-1) + time;
  else
    duration(1) = time;
  end

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

  % Show progress
  fprintf('Step %d: loglike=%e, angle=%e  (%f seconds)\n', k, logP, angle, ...
          time);
  
  % Check whether to save the results
  if (now-lastsave)*3600*24 >= opts.autosavetime && opts.autosavetime > 0
    fprintf('Saving to %s..', opts.autosavefile);
    save(opts.autosavefile, 'W','X','mu','v','tau','u','a_u','b_u','invB', ...
         'duration', 'cost', 'rmse', 'rmse_test');
    lastsave = now;
    fprintf(' done.\n');
  end


  
end

fprintf('Final log-likelihood: %e\n', logP);
loglikelihood = logP;

if nargout == 1
  Wmat = W;
  clear W;
  W.W = Wmat;
  W.X = X;
  W.invB = invB;
  W.mu = mu;
  W.v = v;
  W.tau = tau;
  W.u = u;
  W.a_u = a_u;
  W.b_u = b_u;
  W.loglikelihood = loglikelihood;
end

% $$$ %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% $$$ function logP = loglikelihood_postdebug(Y,W,mu,s2,nu)
% $$$ nsamp = 100;
% $$$ [m,n] = size(Y);
% $$$ usamp = zeros(1,nsamp);
% $$$ logP = 0;
% $$$ for j=1:n
% $$$   pdf = 0;
% $$$   for k=1:nsamp
% $$$     a = 0.5 * (m+nu);
% $$$     b = 0.5 * ((Y(:,j)-mu)' * inv(W*W' + s2.*eye(m)) * (Y(:,j)-mu) + nu);
% $$$     usamp(1,:) = gamrnd(a, 1/b, 1,nsamp);
% $$$     u = usamp(k);
% $$$     V = (W*W'+ eye(m,m).*diag(s2)) / u;
% $$$     pdf = pdf + 1/mnorm_pdf((Y(:,j)-mu)',0,V);
% $$$   end
% $$$   pdf = pdf / nsamp;
% $$$   logP = logP - log(pdf);
% $$$ end
% $$$ 
% $$$ %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% $$$ function logP = loglikelihood_priordebug(Y,W,mu,s2,nu)
% $$$ nsamp = 100;
% $$$ [m,n] = size(Y);
% $$$ usamp = zeros(1,nsamp);
% $$$ logP = 0;
% $$$ for j=1:n
% $$$   usamp(1,:) = gamrnd(nu/2, 2/nu, 1,nsamp);
% $$$   pdf = 0;
% $$$   for k=1:nsamp
% $$$     u = usamp(k);
% $$$     V = (W*W'+ eye(m).*diag(s2)) / u;
% $$$     pdf = pdf + mnorm_pdf((Y(:,j)-mu)',0,V);
% $$$   end
% $$$   pdf = pdf / nsamp;
% $$$   logP = logP + log(pdf);
% $$$ end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [X,CovX,uXX,W,mu] = orthogonalize(X,u,CovX,uXX,W,mu)

m = size(W,1);
[d,n] = size(X);

% "zero mean"
dmu = (X*u') / sum(u);
X = X - repmat(dmu,1,n);
mu = mu + W*dmu;

% "orthogonalize" w.r.t. X
meanUXX = sum(uXX,3) / n;
[Vx,Dx] = eig(meanUXX);
Rx = diag(1./sqrt(diag(Dx))) * Vx';
Rw = Vx * sqrt(Dx);
X = Rx * X;
for j=1:n
  uXX(:,:,j) = Rx * uXX(:,:,j) * Rx';
  CovX(:,:,j) = Rx * CovX(:,:,j) * Rx';
end
W = W * Rw;

% rotate w.r.t. W
[ W, Dw, Vw ] = svd(W);
W = W*Dw;
Qx = Vw';
X = Qx*X;
for j = 1:n
  uXX(:,:,j) = Qx*uXX(:,:,j)*Qx';
  CovX(:,:,j) = Qx * CovX(:,:,j) * Qx';
end
