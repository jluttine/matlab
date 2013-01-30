function [W,CovW,X,CovX] = vbgppcamv(Y,D,inW,inX,covfuncW,logthetaW,covfuncX, ...
                                     logthetaX,varargin)
% [W,CovW,X,CovX] = vbgppcamv(Y,D,inW,inX,covfuncW,logthetaW,covfuncX, ...
% logthetaX,varargin)
% 
% Variational Bayesian (VB) Gaussian process (GP) principal component
% analysis (PCA) with handling of missing values (MV).

[M,N] = size(Y);
% $$$ vecM = 1:M;
% $$$ vecN = 1:N;
% $$$ sizeX = [D,N];
% $$$ sizeW = [M,D];

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
  inW(:,rowrm) = [];
  Y(rowrm,:) = [];
  Obs(rowrm,:) = [];
  Morig = M;
  M = rows(Y);
  W(rowrm,:) = [];
  varW(rowrm,:) = [];
end
if any(colrm)
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

% Number of pseudo inputs for each component
opts.pseudodensityx = opts.pseudodensityx(:) .* ones(D,1);
opts.pseudodensityw = opts.pseudodensityw(:) .* ones(D,1);
Np = ceil(opts.pseudodensityx.*cols(inX)); % no. of pseudos for X components
Mp = ceil(opts.pseudodensityw.*cols(inW)); % no. of pseudos for W components

dimX = rows(inX);   % dimensionality of the input space of X
dimW = rows(inW);   % dimensionality of the input space of W

% Initialize pseudo inputs
pseudoX = cell(D,1);
pseudoW = cell(D,1);
usepseudoX = 1==zeros(D,1);
usepseudoW = 1==zeros(D,1);
for d=1:D
  if Np(d)==N
    pseudoX{d} = inX; % full, don't use pseudo inputs
    usepseudoX(d) = false;
  else
    permN = randperm(N);
    pseudoX{d} = inX(:,permN(1:Np(d)));
    usepseudoX(d) = true;
  end
  if Mp(d)==M
    pseudoW{d} = inW; % full, don't use pseudo inputs
    usepseudoW(d) = false;
  else
    permM = randperm(M);
    pseudoW{d} = inW(:,permM(1:Mp(d)));
    usepseudoW(d) = true;
  end
end

if ~isempty(opts.initpseudow)
  pseudoW = opts.initpseudow;
  usepseudoW(:) = true;
end
if ~isempty(opts.initpseudox)
  pseudoX = opts.initpseudox;
  usepseudoX(:) = true;
end

cholKWpp = cell(D,1);
cholKXpp = cell(D,1);

mu = zeros(M,1);
v_mu = zeros(M,1);

logP = -inf;
costcpu = 0;

% Wp
Wp = cell(D,1);  
CovWp = cell(D,1);

% Xp
Xp = cell(D,1);  
CovXp = cell(D,1);

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

for ind=1:opts.maxiter
  
  oldW = W;
  
  startitercpu = cputime;
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

  % "Whiten" X, a bit adhoc..!! Maybe you should whiten second moment,
  % not variance?!!!
  xx = colsum(X.^2+varX)/N;
  sqrtxx = sqrt(xx);
  %stdX = std(X,1,2);
  varX = bsxfun(@rdivide, varX, xx);
  X = bsxfun(@rdivide, X, sqrtxx);
  varW = bsxfun(@times, varW, xx');
  W = bsxfun(@times, W, sqrtxx');
    
  %% Update W
  for d=1:D
    others = [1:(d-1), (d+1):D];
    V = spalloc(M,M,M);
    c = zeros(M,1);
    for m=1:M
      obs = Obs(m,:);
      V(m,m) = tau * (X(d,obs)*X(d,obs)' + sum(varX(d,obs)));
      c(m) = tau * X(d,obs) * (Y(m,obs) - W(m,others)*X(others,obs))';
    end
    
    if usepseudoW(d)
      pseudos = pseudoW{d};
    else
      pseudos = [];
    end
    [Wp{d}, CovWp{d}, pseudos, logthetaW{d}, tmp, Kwp, cholKWpp{d}] = ...
        gplearn(logthetaW{d}, covfuncW{d}, inW, [], V, pseudos, 'vy', c, ...
                'maxsearch', maxsearchw, 'checkgrad', opts.checkgradw, ...
                'updatepseudo', opts.updatepseudow);

    if ~usepseudoW(d)
      % full
      W(:,d) = Wp{d};
      varW(:,d) = diag(CovWp{d});
    else
      % using pseudo inputs
      pseudoW{d} = pseudos;
      [W(:,d), varW(:,d)] = gppred(pseudoW{d}, Wp{d}, CovWp{d}, inW, ...
                                   logthetaW{d}, covfuncW{d}, 'cholkx', ...
                                   cholKWpp{d}, 'khx', Kwp);
    end
    ltW_in_vbgppca = logthetaW{d};
  end
  

  %% Update X
  for d=1:D
    others = [1:(d-1), (d+1):D];
    V = spalloc(N,N,N);
    c = zeros(N,1);
    for n=1:N
      obs = Obs(:,n);
      V(n,n) = tau * (W(obs,d)'*W(obs,d) + sum(varW(obs,d)));
      c(n) = tau * W(obs,d)' * (Y(obs,n) - W(obs,others)*X(others,n));
    end
    if usepseudoX(d)
      pseudos = pseudoX{d};
    else
      pseudos = [];
    end
    [Xp{d}, CovXp{d}, pseudos, logthetaX{d}, tmp, Kxp, cholKXpp{d}] = ...
        gplearn(logthetaX{d}, covfuncX{d}, inX, [], V, pseudos, 'vy', c, ...
                'maxsearch', maxsearchx, 'checkgrad', opts.checkgradx, ...
                'updatepseudo', opts.updatepseudox);

    if ~usepseudoX(d)
      % full
      X(d,:) = Xp{d};
      varX(d,:) = diag(CovXp{d});
    else
      % pseudo inputs
      pseudoX{d} = pseudos;
      [X(d,:), varX(d,:)] = gppred(pseudoX{d}, Xp{d}, CovXp{d}, inX, ...
                                   logthetaX{d}, covfuncX{d}, 'cholkx', ...
                                   cholKXpp{d}, 'khx', Kxp);
    end
    
    ltX_in_vbgppca = logthetaX{d};
  end
  
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
           + traceprod( (W(obs,:)'*W(obs,:) + diag(rowsum(varW(obs,:))) ), ...
                        (X(:,n)*X(:,n)' + diag(varX(:,n))) ) ...
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
  
    % Cost from Y
    cost_Y = -0.5*tau*err2 - 0.5*NMobs*(log2pi - logtau);

    % Cost from W
    cost_W = 0;
    for d=1:D
      cost_W = cost_W + mvnvbcost(Wp{d},CovWp{d}, 0,0, cholKWpp{d}, ...
                                  2*logdettri(cholKWpp{d}));
    end

    % Cost from X
    cost_X = 0;
    for d=1:D
      cost_X = cost_X + mvnvbcost(Xp{d},CovXp{d}, 0,0, cholKXpp{d}, ...
                                  2*logdettri(cholKXpp{d}));
    end

    % Cost from tau
    cost_tau = -gamkl(tau,logtau,a_tau,b_tau,atau,btau);
    
    oldlogP = logP;

    % Lower bound of loglikelihood
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


% Predict whole dataset, also rows and columns that were removed at the
% beginning because of they were empty
if opts.reconstruct
  if any(rowrm)
    W = zeros(Morig,D);
    varW = zeros(Morig,D);
    for d=1:D
      [W(:,d), varW(:,d)] = gppred(pseudoW{d}, Wp{d}, CovWp{d}, inWorig, ...
                                   logthetaW{d}, covfuncW{d});
    end
  end
  if any(colrm)
    X = zeros(D,Norig);
    varX = zeros(D,Norig);
    for d=1:D
      [X(d,:), varX(d,:)] = gppred(pseudoX{d}, Xp{d}, CovXp{d}, inXorig, ...
                                   logthetaX{d}, covfuncX{d});
    end    
  end
end

if nargout == 1
  Q.Wp = Wp;
  Q.CovWp = CovWp;
  Q.W = W;
  Q.varW = varW;
  Q.pseudoW = pseudoW;
  Q.logthetaW = logthetaW;
  Q.covfuncW = covfuncW;
  Q.inW = inWorig;
  
  Q.Xp = Xp;
  Q.CovXp = CovXp;
  Q.X = X;
  Q.varX = varX;
  Q.pseudoX = pseudoX;
  Q.logthetaX = logthetaX;
  Q.covfuncX = covfuncX;
  Q.inX = inXorig;
  
  Q.tau = tau;
  Q.a_tau = a_tau;
  Q.b_tau = b_tau;
  
  Q.loglikelihood = loglikelihood;
  
  W = Q;
end


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
  varW = zeros(M,D);
end

% X
if isstruct(init) && isfield(init, 'X') && ~isempty(init.X)
  disp('Initialize X to given value');
  X = init.X;
else
  %  X = randn(D,N); 
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
  varX = 0 * ones(D,N);
end

%% tau precision (inverse noise)
if isstruct(init) && isfield(init, 'tau') && ~isempty(init.tau)
  tau = init.tau;
else
  tau = 1e3;
end
