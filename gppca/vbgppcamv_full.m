function [W,CovW,X,CovX] = vbgppcamv_full(Y,D,inW,inX,covfuncW,logthetaW, ...
                                          covfuncX, logthetaX,varargin)
% [W,CovW,X,CovX] = vbgppcamv(Y,D,inW,inX,covfuncW,logthetaW,covfuncX, ...
% logthetaX,varargin)
% 
% Variational Bayesian (VB) Gaussian process (GP) principal component
% analysis (PCA) with handling of missing values (MV).
%
% This version does NOT factorize with respect to the components.
%
% At the moment, pseudo inputs are not yet supported.
%

[M,N] = size(Y);

opts = struct( ...
    'init',          [],...
    'rotate', true, ...
    'autosavetime', 0,...
    'autosavefile', 'vbgppcamv_autosave',...
    'testset', [], ...
    'pseudodensityx', 1, ...
    'pseudodensityw', 1, ...
    'updatepseudox', true, ...
    'updatepseudow', true, ...
    'initpseudow', [], ...
    'initpseudox', [], ...
    'reconstruct', true, ...
    'loglikelihood', true, ...
    'updatehyper', 1, ...
    'maxsearchx', 10, ...
    'maxsearchw', 10, ...
    'checkgradx', false, ...
    'checkgradw', false, ...
    'maxiter', 100);

[ opts, errmsg, wrnmsg ] = argschk( opts, varargin{:} );
if ~isempty(errmsg), error( errmsg ), end
if ~isempty(wrnmsg), warning( wrnmsg ), end

Obs = ~isnan(Y);
NMobs = sum(Obs(:));

inWorig = inW;
inXorig = inX;

% Initialize posterior values
[W,varW,X,varX,tau] = initialize(opts.init, D, inW,inX,covfuncW, ...
                                 logthetaW,covfuncX, logthetaX);

% Remove empty rows/columns
rowrm = (colsum(Obs) == 0);
colrm = (rowsum(Obs) == 0);
if any(rowrm)
  disp('Removing empty rows');
  inW(:,rowrm) = [];
  Y(rowrm,:) = [];
  Obs(rowrm,:) = [];
  Morig = M;
  M = rows(Y);
  W(rowrm,:) = [];
  varW(rowrm,:) = [];
end
if any(colrm)
  disp('Removing empty columns');
  inX(:,colrm) = [];
  Y(:,colrm) = [];
  Obs(:,colrm) = [];
  Norig = N;
  N = cols(Y);
  X(:,colrm) = [];
  varX(:,colrm) = [];
end

log2pi = log(2*pi);


%%
%% Initialize posterior parameters

% $$$ % Number of pseudo inputs for each component
% $$$ opts.pseudodensityx = opts.pseudodensityx(:) .* ones(D,1);
% $$$ opts.pseudodensityw = opts.pseudodensityw(:) .* ones(D,1);
% $$$ Np = ceil(opts.pseudodensityx.*cols(inX)); % no. of pseudos for X components
% $$$ Mp = ceil(opts.pseudodensityw.*cols(inW)); % no. of pseudos for W components
% $$$ 
% $$$ dimX = rows(inX);   % dimensionality of the input space of X
% $$$ dimW = rows(inW);   % dimensionality of the input space of W
% $$$ 
% $$$ % Initialize pseudo inputs
% $$$ pseudoX = cell(D,1);
% $$$ pseudoW = cell(D,1);
% $$$ usepseudoX = 1==zeros(D,1);
% $$$ usepseudoW = 1==zeros(D,1);
% $$$ for d=1:D
% $$$   if Np(d)==N
% $$$     pseudoX{d} = inX; % full, don't use pseudo inputs
% $$$     usepseudoX(d) = false;
% $$$   else
% $$$     permN = randperm(N);
% $$$     pseudoX{d} = inX(:,permN(1:Np(d)));
% $$$     usepseudoX(d) = true;
% $$$   end
% $$$   if Mp(d)==M
% $$$     pseudoW{d} = inW; % full, don't use pseudo inputs
% $$$     usepseudoW(d) = false;
% $$$   else
% $$$     permM = randperm(M);
% $$$     pseudoW{d} = inW(:,permM(1:Mp(d)));
% $$$     usepseudoW(d) = true;
% $$$   end
% $$$ end
% $$$ 
% $$$ if ~isempty(opts.initpseudow)
% $$$   pseudoW = opts.initpseudow;
% $$$   usepseudoW(:) = true;
% $$$ end
% $$$ if ~isempty(opts.initpseudox)
% $$$   pseudoX = opts.initpseudox;
% $$$   usepseudoX(:) = true;
% $$$ end
% $$$ 
% $$$ cholKWpp = cell(D,1);
% $$$ cholKXpp = cell(D,1);

mu = zeros(M,1);
v_mu = zeros(M,1);

logP = -inf;
costcpu = 0;

% $$$ % Wp
% $$$ Wp = cell(D,1);  
% $$$ CovWp = cell(D,1);
% $$$ 
% $$$ % Xp
% $$$ Xp = cell(D,1);  
% $$$ CovXp = cell(D,1);

%% tau precision (inverse noise)
atau = 1e-4;
btau = 1e-4;
a_tau = nan; % posterior pdf parameter
b_tau = nan; % posterior pdf parameter
logtau = nan;

lastsave = cputime;

%%
%% VB learning
% $$$ K = myfeval(covfuncW{1}{:}, logthetaW{1}, inW, inW);
% $$$ figure
% $$$ imagesc(K)
% $$$ return

% Index order (component-wise, i.e., block-diagonal K)
indsW = reshape(1:(D*M), [M D]);
indsX = reshape(1:(D*N), [N D])';

CovW = eye(D*M);
CovX = eye(D*N);

% $$$ %% TEMP
% $$$ W(:) = 1;
% $$$ CovW(:) = 0;
 
% Initialize help variables
WW = zeros(D,D,M);
for m=1:M
  WW(:,:,m) = W(m,:)'*W(m,:) + CovW(indsW(m,:),indsW(m,:));
end
XX = zeros(D,D,N);
for n=1:N
  XX(:,:,n) = X(:,n)*X(:,n)' + CovX(indsX(:,n),indsX(:,n));
end

% Evaluate prior covariance matrices
KW = zeros(D*M);
KX = zeros(D*N);
for d=1:D
  first = (d-1)*M + 1;
  last = first + M - 1;
  inds = first:last;
  covfunc = covfuncW{d};
  if ~iscell(covfunc)
    covfunc = {covfunc};
  end
  KW(inds,inds) = feval(covfunc{:}, logthetaW{d}, inW, inW);
  
  first = (d-1)*N + 1;
  last = first + N - 1;
  inds = first:last;
  covfunc = covfuncX{d};
  if ~iscell(covfunc)
    covfunc = {covfunc};
  end
  KX(inds,inds) = feval(covfunc{:}, logthetaX{d}, inX, inX);
end
% $$$ KW = regularize(KW, 1e1);
% $$$ LKW = chol(KW, 'lower');
% $$$ KX = regularize(KX, 1e2);
% $$$ LKX = chol(KX, 'lower');

% $$$ figure
% $$$ imagesc(KW)
% $$$ figure
% $$$ imagesc(KX)
  

reg = 1e-6;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% VB EM

for ind=1:opts.maxiter
  
  oldW = W;
  
  startitercpu = cputime;
  
  %% UPDATE VARIABLES
  
  update = false;
  if length(opts.updatehyper) > 1 
    if any(ind==opts.updatehyper)
      update = true;
    end
  else
    if mod(ind, opts.updatehyper) == 0 && opts.updatehyper ~=0
      update = true;
    end
  end
  if update
    maxsearchx = opts.maxsearchx;
    maxsearchw = opts.maxsearchw;
  else
    maxsearchx = 0;
    maxsearchw = 0;
  end

  if ind < 5 && false
    % Scale X to unit variance. I think it should improve the lower bound
    % but it doesn't seem to do so in practice (at least always).. ??
    disp('Using adhoc whitening of second moment')
    %xx = diag(sum(XX,3))/N - mean(X,2).^2;   % variance
    xx = diag(sum(XX,3))/N;   % second moment
    sqrtxx = sqrt(xx);
    X = bsxfun(@rdivide, X, sqrtxx);
    XX = bsxfun(@rdivide, XX, (sqrtxx(:)*sqrtxx(:)'));
    W = bsxfun(@times, W, sqrtxx');
    WW = bsxfun(@times, WW, (sqrtxx(:)*sqrtxx(:)'));
  end
  %XX_norm = mean(XX,3)
    
  %% Update W
  U = spalloc(D*M,D*M, M*D*D);
  z = zeros(D*M,1);
  for m=1:M
    obs = Obs(m,:);
    U(indsW(m,:),indsW(m,:)) = tau * sum(XX(:,:,obs),3);
    z(indsW(m,:)) = tau * X(:,obs) * Y(m,obs)';
  end
  % Optimise hyperparameters
  if update
    logthetaW = optimise_hyperparameter(logthetaW, covfuncW, inW, U, z, D, M);
  end
  % Update prior covariance matrix
  KW = zeros(D*M);
  for d=1:D
    first = (d-1)*M + 1;
    last = first + M - 1;
    inds = first:last;
    covfunc = covfuncW{d};
    if ~iscell(covfunc)
      covfunc = {covfunc};
    end
    KW(inds,inds) = feval(covfunc{:}, logthetaW{d}, inW, inW);
    
    KW = KW + reg*eye(size(KW));
    
  end
  %KW = regularize(KW, 1e1);
  %LKW = chol(KW, 'lower');
  % Evaluate posterior
  if false % choose method
    L = chol(KW + KW*U*KW, 'lower');
    sqrtCovW = solve_tril(L, KW);
    CovW = sqrtCovW' * sqrtCovW; 
    W(:) = CovW(indsW(:),indsW(:)) * z(indsW(:));
  else
    % Prefer this!
%    L = safechol(KW + inv(U), 1e-8, 'lower');
    L = chol(KW + inv(U), 'lower');
    tmp = solve_tril(L, KW);
    CovW = KW - tmp' * tmp;
    tmp = linsolve_chol(L, U \ z, 'lower');
    W(:) = KW(indsW(:),indsW(:)) * tmp(indsW(:));
  end
  WW = zeros(D,D,M);
  for m=1:M
    WW(:,:,m) = W(m,:)'*W(m,:) + CovW(indsW(m,:),indsW(m,:));
  end
% $$$   figure
% $$$   imagesc(WW(:,:,1))
% $$$   error('jou')
  
% $$$   figure
% $$$   imagesc(CovW)
% $$$   error('halo')
  
  %% Update X
  U = spalloc(D*N,D*N, N*D*D);
  z = zeros(D*N,1);
  for n=1:N
    obs = Obs(:,n);
    U(indsX(:,n),indsX(:,n)) = tau * sum(WW(:,:,obs),3);
    z(indsX(:,n)) = tau * W(obs,:)' * Y(obs,n);
  end
  % Optimise hyperparameters
  if update
    logthetaX = optimise_hyperparameter(logthetaX, covfuncX, inX, U, z, D, N);
  end
  % Update prior covariance matrix
  KX = zeros(D*N);
  for d=1:D
    first = (d-1)*N + 1;
    last = first + N - 1;
    inds = first:last;
    covfunc = covfuncX{d};
    if ~iscell(covfunc)
      covfunc = {covfunc};
    end
    KX(inds,inds) = feval(covfunc{:}, logthetaX{d}, inX, inX);
  
    KX = KX + reg*eye(size(KX));
      
  end
  %KX = regularize(KX, 1e2);
  %LKX = chol(KX, 'lower');
  % Evaluate posterior
% $$$   figure
% $$$   imagesc(inv(U))
% $$$   figure
% $$$   imagesc(KX)
% $$$   tmp = chol(U, 'lower');
% $$$   figure
% $$$   imagesc(linsolve_chol(tmp, eye(size(tmp))) == 0);
  if false
    L = chol(KX + KX*U*KX, 'lower');
    sqrtCovX = solve_tril(L, KX);
    CovX = sqrtCovX' * sqrtCovX; 
    X(:) = CovX(indsX(:),indsX(:)) * z(indsX(:));
  else
    % Prefer this!
    L = chol(KX + inv(U), 'lower');
%    L = safechol(KX + inv(U), 1e-8, 'lower');
    tmp = solve_tril(L, KX);
    CovX = KX - tmp' * tmp;
    tmp = linsolve_chol(L, U \ z, 'lower');
    X(:) = KX(indsX(:),indsX(:)) * tmp(indsX(:));
  end
  XX = zeros(D,D,N);
  for n=1:N
    XX(:,:,n) = X(:,n)*X(:,n)' + CovX(indsX(:,n),indsX(:,n));
  end
  
  %XX_norm = mean(XX,3)

%  CovX(1:10,1:10)
%  tmp = chol(KX, 'lower');

% $$$   figure
% $$$   imagesc(CovX)
% $$$   error('halo')
  
% $$$   rcondX = rcond(CovX)
% $$$   rcondestX = rcond(CovX + 1e-4*normest(CovX)*eye(size(CovX)))
% $$$   tmp = chol(CovX,'lower');
  
%  [X,CovX,W,CovW] = rotate(X,CovX,indsX,LKX,W,CovW);
  
% $$$   for d=1:D
% $$$     others = [1:(d-1), (d+1):D];
% $$$     V = spalloc(M,M,M);
% $$$     c = zeros(M,1);
% $$$     for m=1:M
% $$$       obs = Obs(m,:);
% $$$       V(m,m) = tau * (X(d,obs)*X(d,obs)' + sum(varX(d,obs)));
% $$$       c(m) = tau * X(d,obs) * (Y(m,obs) - W(m,others)*X(others,obs))';
% $$$     end
% $$$     
% $$$     if usepseudoW(d)
% $$$       pseudos = pseudoW{d};
% $$$     else
% $$$       pseudos = [];
% $$$     end
% $$$     [Wp{d}, CovWp{d}, pseudos, logthetaW{d}, tmp, Kwp, cholKWpp{d}] = ...
% $$$         gplearn(logthetaW{d}, covfuncW{d}, inW, [], V, pseudos, 'vy', c, ...
% $$$                 'maxsearch', maxsearchw, 'checkgrad', opts.checkgradw, ...
% $$$                 'updatepseudo', opts.updatepseudow);
% $$$ 
% $$$     if ~usepseudoW(d)
% $$$       % full
% $$$       W(:,d) = Wp{d};
% $$$       varW(:,d) = diag(CovWp{d});
% $$$     else
% $$$       % using pseudo inputs
% $$$       pseudoW{d} = pseudos;
% $$$       [W(:,d), varW(:,d)] = gppred(pseudoW{d}, Wp{d}, CovWp{d}, inW, ...
% $$$                                    logthetaW{d}, covfuncW{d}, 'cholkx', ...
% $$$                                    cholKWpp{d}, 'khx', Kwp);
% $$$     end
% $$$     ltW_in_vbgppca = logthetaW{d};
% $$$   end
  

% $$$   %% Update X
% $$$   for d=1:D
% $$$     others = [1:(d-1), (d+1):D];
% $$$     V = spalloc(N,N,N);
% $$$     c = zeros(N,1);
% $$$     for n=1:N
% $$$       obs = Obs(:,n);
% $$$       V(n,n) = tau * (W(obs,d)'*W(obs,d) + sum(varW(obs,d)));
% $$$       c(n) = tau * W(obs,d)' * (Y(obs,n) - W(obs,others)*X(others,n));
% $$$     end
% $$$     if usepseudoX(d)
% $$$       pseudos = pseudoX{d};
% $$$     else
% $$$       pseudos = [];
% $$$     end
% $$$     [Xp{d}, CovXp{d}, pseudos, logthetaX{d}, tmp, Kxp, cholKXpp{d}] = ...
% $$$         gplearn(logthetaX{d}, covfuncX{d}, inX, [], V, pseudos, 'vy', c, ...
% $$$                 'maxsearch', maxsearchx, 'checkgrad', opts.checkgradx, ...
% $$$                 'updatepseudo', opts.updatepseudox);
% $$$ 
% $$$     if ~usepseudoX(d)
% $$$       % full
% $$$       X(d,:) = Xp{d};
% $$$       varX(d,:) = diag(CovXp{d});
% $$$     else
% $$$       % pseudo inputs
% $$$       pseudoX{d} = pseudos;
% $$$       [X(d,:), varX(d,:)] = gppred(pseudoX{d}, Xp{d}, CovXp{d}, inX, ...
% $$$                                    logthetaX{d}, covfuncX{d}, 'cholkx', ...
% $$$                                    cholKXpp{d}, 'khx', Kxp);
% $$$     end
% $$$     
% $$$     ltX_in_vbgppca = logthetaX{d};
% $$$   end
  
  disp('Updating tau')
  
  %% Update tau
  % TODO: You could loop over the longer dimension N or M? Then
  % vectorization would help the most. :)
  err2 = 0;
  for n=1:N
    obs = Obs(:,n);
    ymu = Y(obs,n) - mu(obs);
    err2 = err2 ...
           + ymu'*ymu ...
           + traceprod(sum(WW(:,:,obs),3), XX(:,:,n)) ...
           + sum(v_mu) ...
           - 2*(W(obs,:)*X(:,n))' * ymu;
  end
  a_tau = atau + 0.5*NMobs;
  b_tau = btau + 0.5*err2;
  tau = a_tau / b_tau;
  logtau = psi(a_tau) - log(b_tau);
  
  % CPU time used for calculations
  itercpu = cputime - startitercpu;
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% CALCULATE LOWER BOUND OF THE LOG-LIKELIHOOD
  
  disp('Evaluating lower bound.');

  if opts.loglikelihood

    startcostcpu = cputime;
    oldlogP = logP;
  
    % Cost from Y
    cost_Y = -0.5*tau*err2 - 0.5*NMobs*(log2pi - logtau);
    
    % WTF!?!?!?!!! If using safechol, the coefficient has a HUUUGE effect
    % on the loglikelihood lower bound!!!

    % Cost from W
    %rcond_KW = rcond(KW)
%    LKW = chol(KW(indsW(:),indsW(:)), 'lower');
    LKW = chol(KW(indsW(:),indsW(:))+reg*eye(size(KW)), 'lower');
%    LKW = safechol(KW(indsW(:),indsW(:)), 1e-15, 'lower');
    cost_W = mvnvbcost(W(:),CovW(indsW(:),indsW(:)), 0,0, LKW, 2*logdettri(LKW));
             
% $$$     cost_W = 0;
% $$$     for d=1:D
% $$$       cost_W = cost_W + mvnvbcost(Wp{d},CovWp{d}, 0,0, cholKWpp{d}, ...
% $$$                                   2*logdettri(cholKWpp{d}));
% $$$     end

    % Cost from X
%    rcond_KX = rcond(KX)
%    LKX = chol(KX(indsX(:),indsX(:)), 'lower');
    LKX = chol(KX(indsX(:),indsX(:))+reg*eye(size(KX)), 'lower');
%    LKX = safechol(KX(indsX(:),indsX(:)), 1e-15, 'lower');
    cost_X = mvnvbcost(X(:),CovX(indsX(:),indsX(:)), 0,0, LKX, 2*logdettri(LKX));
% $$$     cost_X = 0;
% $$$     for d=1:D
% $$$       cost_X = cost_X + mvnvbcost(Xp{d},CovXp{d}, 0,0, cholKXpp{d}, ...
% $$$                                   2*logdettri(cholKXpp{d}));
% $$$     end

    % Cost from tau
    cost_tau = -gamkl(tau,logtau,a_tau,b_tau,atau,btau);
% $$$     
% $$$     oldlogP = logP;
% $$$ 
% $$$     % Lower bound of loglikelihood
    logP = cost_Y + cost_W + cost_X + cost_tau;% + cost_mu + cost_tau + cost_w;
    costcpu = cputime - startcostcpu;

    if logP < oldlogP
      warning('Cost increased!!');
    end
  
  end

  
  disp('Monitoring stuff.')
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% MONITORING STUFF
  
  % Change in subspace
  angle = subspace(oldW, W);
  
  %% Print progress
  fprintf('%d step: loglike=%e, angle=%.2e (itercpu=%.2fs, nllcpu=%.2fs)\n', ...
          ind, logP, angle, itercpu, costcpu);
  
  %% Save
  loglikelihood = logP;
  if opts.autosavetime > 0 && cputime - lastsave > opts.autosavetime
    fprintf('Saving to a file %s..', opts.autosavefile);
    iteration = ind;
    save(opts.autosavefile, 'Wp', 'CovWp', 'W', 'varW', 'pseudoW', ...
         'logthetaW', 'covfuncW', 'inW', 'Xp', 'CovXp', 'X', 'varX', ...
         'pseudoX', 'logthetaX', 'covfuncX', 'inX', 'tau', 'a_tau', 'b_tau', ...
         'loglikelihood', 'ind');
    fprintf(' done\n');
    lastsave = cputime;
  end
  
end

% $$$ figure
% $$$ imagesc(CovX == 0)


% Predict whole dataset, also rows and columns that were removed at the
% beginning because of they were empty
if opts.reconstruct
  if any(rowrm)
    % Change variables
    W_old = W;
    KW_old = KW;
    CovW_old = CovW;
    Morig = size(inWorig,2);
    W = zeros(Morig,D);
    KW = zeros(D*Morig);
    KWWold = zeros(D*Morig,D*M);
    CovW = zeros(D*Morig);
    indsWorig = reshape(1:(D*Morig), [Morig D]);
    % New prior covariance matrix
    for d=1:D
      first = (d-1)*Morig + 1;
      last = first + Morig - 1;
      inds = first:last;
      first = (d-1)*Morig + 1;
      last = first + Morig - 1;
      inds_old = first:last;
      covfunc = covfuncW{d};
      if ~iscell(covfunc)
        covfunc = {covfunc};
      end
      KWWold(inds,inds_old) = feval(covfunc{:}, logthetaW{d}, inWorig, inW);
      KW(inds,inds) = feval(covfunc{:}, logthetaW{d}, inWorig, inWorig);
      KW = KW + reg*eye(size(KW));
    end
    % New posterior quantities
    LKW_old = chol(KW_old(indsW(:),indsW(:)), 'lower');
    W(:) = KWWold(indsWorig(:),indsW(:)) * linsolve_chol(LKW_old, W_old(:), ...
                                                      'lower');
    A = linsolve_chol(LKW_old, KWWold(indsWorig(:),indsW(:))', 'lower');
    CovW(indsWorig(:),indsWorig(:)) = KW(indsWorig(:),indsWorig(:)) - ...
        A'*(KW_old(indsW(:),indsW(:))-CovW(indsW(:),indsW(:)))*A;

    M = Morig;
    indsW = indsWorig;
    inW = inWorig;
  end
  if any(colrm)
% $$$     figure
% $$$     imagesc(CovX)
    
    % Change variables
    X_old = X;
    KX_old = KX;
    CovX_old = CovX;
    Norig = size(inXorig,2);
    X = zeros(D,Norig);
    KX = zeros(D*Norig);
    KXXold = zeros(D*Norig,D*N);
    CovN = zeros(D*Norig);
    indsXorig = reshape(1:(D*Norig), [Norig D])';
    % New prior covariance matrix
    for d=1:D
% $$$       first = (d-1)*Morig + 1;
% $$$       last = first + Morig - 1;
% $$$       inds = first:last;
      first = (d-1)*Norig + 1;
      last = first + Norig - 1;
      inds = first:last;
      first = (d-1)*N + 1;
      last = first + N - 1;
      inds_old = first:last;
      covfunc = covfuncX{d};
      if ~iscell(covfunc)
        covfunc = {covfunc};
      end
      KXXold(inds,inds_old) = feval(covfunc{:}, logthetaX{d}, inXorig, inX);
      KX(inds,inds) = feval(covfunc{:}, logthetaX{d}, inXorig, inXorig);
      KX = KX + reg*eye(size(KX));
    end
    % New posterior quantities
    LKX_old = chol(KX_old(indsX(:),indsX(:)), 'lower');
    X(:) = KXXold(indsXorig(:),indsX(:)) * linsolve_chol(LKX_old, X_old(:), ...
                                                      'lower');
    A = linsolve_chol(LKX_old, KXXold(indsXorig(:),indsX(:))', 'lower');
    CovX(indsXorig(:),indsXorig(:)) = KX(indsXorig(:),indsXorig(:)) - ...
        A'*(KX_old(indsX(:),indsX(:))-CovX(indsX(:),indsX(:)))*A;
  
% $$$     figure
% $$$     imagesc(CovX)
    
    N = Norig;
    indsX = indsXorig;
    inX = inXorig;
  end
end
% $$$ if opts.reconstruct
% $$$   if any(rowrm)
% $$$     oldW = W;
% $$$     KW_old = KW;
% $$$     KW = feval(covfunc{:}, logthetaW{d}, inWorig, inWorig);
% $$$     W = zeros(Morig,D);
% $$$     varW = zeros(Morig,D);
% $$$     for d=1:D
% $$$       [W(:,d), varW(:,d)] = gppred(inW, oldW(:,d), CovW(indsW(:,d),indsW(:,d)), ...
% $$$                                    inWorig, logthetaW{d}, covfuncW{d});
% $$$     end
% $$$   end
% $$$   if any(colrm)
% $$$     oldX = X;
% $$$     X = zeros(D,Norig);
% $$$     varX = zeros(D,Norig);
% $$$     for d=1:D
% $$$       [X(d,:), varX(d,:)] = gppred(inX, oldX(d,:)', CovX(indsX(d,:),indsX(d,:)), ...
% $$$                                    inXorig, logthetaX{d}, covfuncX{d});
% $$$     end    
% $$$   end
% $$$ end

varX = reshape(diag(CovX), [N D])';
varW = reshape(diag(CovW), [M D]);

if nargout == 1
  Q.W = W;
  Q.CovW = CovW;
% $$$   Q.W = W;
  Q.varW = varW;
% $$$   Q.pseudoW = pseudoW;
  Q.logthetaW = logthetaW;
  Q.covfuncW = covfuncW;
  Q.inW = inWorig;
  
  Q.X = X;
  Q.CovX = CovX;
% $$$   Q.X = X;
  Q.varX = varX;
% $$$   Q.pseudoX = pseudoX;
  Q.logthetaX = logthetaX;
  Q.covfuncX = covfuncX;
  Q.inX = inXorig;
  
  Q.tau = tau;
  Q.a_tau = a_tau;
  Q.b_tau = b_tau;
  
  Q.indsW = indsW;
  Q.indsX = indsX;
  
  Q.loglikelihood = loglikelihood;
  
  W = Q;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function logtheta = optimise_hyperparameter(logtheta, covfuncs, input, U, ...
                                            z, D, N)

% Because the cell logtheta must be presented as one long vector, mark
% the indeces for each component. That is, the hyperparameters of d-th
% component are vector( (indeces(d-1)+1):indeces(d) ) = logtheta{d}
vector = [];
indeces = zeros(D,1);
for d=1:D
  vector = [vector; logtheta{d}(:)];
  indeces(d) = length(logtheta{d});
end
indeces = cumsum(indeces);

invU = inv(U);

%hyperparameters_before = exp(vector(:))

% $$$ mycheckgrad(@hyperparameter_upperbound, vector, 1e-9, covfuncs, input, ...
% $$$             invU, z, D, N, indeces);

vector = minimize(vector, @hyperparameter_upperbound, 5, covfuncs, input, ...
                  invU, z, D, N, indeces);

%hyperparameters_after = exp(vector(:))

ind = 1;
for d=1:D
  logtheta{d} = vector(ind:indeces(d));
  ind = indeces(d) + 1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [f, df] = hyperparameter_upperbound(logtheta_vector, covfuncs, ...
                                             input, invU, z, D, N, ...
                                             logtheta_indeces)

% Projected pseudo-observations
y = invU * z;

% Transform the hyperparameter vector to cell
ind = 1;
for d=1:D
  logtheta{d} = logtheta_vector(ind:logtheta_indeces(d));
  ind = logtheta_indeces(d) + 1;
end

dlogtheta = cell(D,1);

K = zeros(D*N);
for d=1:D
  % Indeces of a component block
  first = (d-1)*N + 1;
  last = first + N - 1;
  inds = first:last;
  % Covariance matrix for the block
  covfunc = covfuncs{d};
  if ~iscell(covfunc)
    covfunc = {covfunc};
  end
  [K(inds,inds), dlogtheta{d}] = feval(covfunc{:}, logtheta{d}, input, input);
end
%K = regularize(K, 1e2);
%LK = chol(K, 'lower');
%K = 

%if rcond(K) < 1e-12
%  f = -

L = chol(K + invU, 'lower');

% Loglikelihood lower bound (mnorm_lpdf)
%f = mnorm_lpdf((invU*z)', 0, K + invU);
v = solve_tril(L, (y-0));
f = -0.5*length(L)*log(2*pi) - logdettri(L) - 0.5 * v'*v;

% Gradients
df = [];
% TODO: USE LINSOLVE TO OPTIMISE BACKSUBSTITUTIONS!!!!
%invLL = L' \ (L \ eye(size(L)));
%c = invLL * y;
c = linsolve_chol(L, y, 'lower');
for d=1:D
  % Indeces of a component block
  first = (d-1)*N + 1;
  last = first + N - 1;
  inds = first:last;
  Ld = L(inds,inds);
  invL = linsolve_chol(Ld, eye(size(Ld)), 'lower');
  % Evaluate gradients
  for i=1:length(logtheta{d})
    dK = dlogtheta{d}(:,:,i);
    %sz_c = size(c)
    %sz_dK = size(dK)
    df = [df; 0.5*c(inds)'*dK*c(inds) - 0.5*traceprod(invL,dK)];
  end
end

% Lower bound to upper bound (because using minimize)
f = -f;
df = -df;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [X,CovX,W,CovW] = rotate(X,CovX,indsX,LKX,W,CovW)

error('Ei taa oikein toimi, ei kannata kayttaa..')

disp('Finding rotation..')
[D,N] = size(X);
R = orth(randn(D));
R = R*R';

%mycheckgrad(@rotation_cost, R(:), 1e-6, X,CovX,indsX,LKX);


R = minimize(R(:), @rotation_cost, 3, X,CovX,indsX,LKX);

R = reshape(R,[D,D]);

% $$$ % Force positive definiteness
% $$$ [VR,DR] = eig(R);
% $$$ R = VR * exp(DR) * inv(VR);

X = R * X;
kronR = kron(R,eye(N));
CovX = kronR * CovX * kronR';

M = rows(W);
invR = inv(R);
W = W * invR;
kronR = kron(invR,eye(M));
CovW = kronR' * CovW * kronR;
disp('Rotation done.');

R
%error('done')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [f,df] = rotation_cost(R,X,CovX,indsX,LKX)
[D,N] = size(X);
R = reshape(R,[D D]);

% Force positive definiteness
% $$$ [VR,DR] = eig(R);
% $$$ R = VR * exp(DR) * inv(VR);

log_qX = -N * logdet(R);
if ~isreal(log_qX)
  f = inf;
  df = nan*ones(D*D,1);
  return
end

LinvKX = solve_tril(LKX, eye(D*N));
invKX = LinvKX' * LinvKX;
%invKX = solve_triu(LKX', solve_tril(LKX, eye(D*N)));
H = 0;
for n=1:N
  for m=1:N
    in = indsX(:,n);
    im = indsX(:,m);
    H = H + invKX(in,im) * R * (X(:,m)*X(:,n)' + CovX(im,in));
  end
end
error('Check that traceprod')
log_pX = -0.5 * traceprod(H,R);

f = -(log_pX - log_qX)
df = N*inv(R') - H;
df = -df(:);

% Test:
X = R * X;
kronR = kron(R,eye(N));
CovX = kronR * CovX * kronR';
cost_X = mvnvbcost(X(:),CovX(indsX(:),indsX(:)), 0,0, LKX, 2*logdettri(LKX))


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ind = row2ind(M,N,m)
ind = (0:(N-1)) * M + m;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ind = col2ind(M,N,n)
ind = (1:M) + (n-1)*M;

%%%%%%%%%%%%%%%%%%%%%
function Y = gpinv(X)
[Q,R] = qr(X);
Y = inv(R) * Q';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [W,varW,X,varX,tau] = initialize(init,D,inW,inX,covfuncW, ...
                                          logthetaW,covfuncX, logthetaX)

M = cols(inW);
N = cols(inX);

% W
if isstruct(init) && isfield(init, 'W') && ~isempty(init.W)
  disp('Initialize W to given value');
  W = init.W;
else
%  W = randn(M,D); 
  disp('Initialize W to zero');%by sampling');
  W = zeros(M,D);
% $$$   for d=1:D
% $$$     W(:,d) = gprnd(inW, logthetaW{d}, covfuncW{d});
% $$$   end
end

% varW
if isstruct(init) && isfield(init, 'varW') && ~isempty(init.varW)
  varW = init.varW;
else
  varW = 1*ones(M,D);
end

% X
if isstruct(init) && isfield(init, 'X') && ~isempty(init.X)
  disp('Initialize X to given value');
  X = init.X;
else
  %X = randn(D,N); 
  disp('Initialize X by sampling');
  X = zeros(D,N);
  for d=1:D
    X(d,:) = gprnd(inX, logthetaX{d}, covfuncX{d});
  end
end

% varX
if isstruct(init) && isfield(init, 'varX') && ~isempty(init.varX)
  varX = init.varX;
else
  varX = 1 * ones(D,N);
end

%% tau precision (inverse noise)
if isstruct(init) && isfield(init, 'tau') && ~isempty(init.tau)
  tau = init.tau;
else
  tau = 1e3;
end
