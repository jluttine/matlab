function [W,CovW,X,CovX] = vbgppcamv(Y,D,inW,inX,kW,thetaW,kX,thetaX,varargin)
% Variational Bayesian (VB) Gaussian process (GP) principal component
% analysis (PCA) with handling of missing values (MV).


[M,N] = size(Y);
vecM = 1:M;
vecN = 1:N;
% $$$ sizeX = [D,N];
% $$$ sizeW = [M,D];

opts = struct( ...
    'init',          [],...
    'rotate', true, ...
    'autosavetime', 0,...
    'autosavefile', 'vbgppcamv_autosave',...
    'testset', [], ...
    'maxiter', 100, ...
    'inputW', 1:M, ...
    'inputX', 1:N, ...
    'covfunX', @(in1,in2,theta) gpcov_ratquad(gpdist(in1,in2),theta(1), ...
                                              theta(2),theta(3)));

[ opts, errmsg, wrnmsg ] = argschk( opts, varargin{:} );
if ~isempty(errmsg), error( errmsg ), end
if ~isempty(wrnmsg), warning( wrnmsg ), end


Obs = ~isnan(Y);

% Form the GP prior covariance matrices using covariance functions
Kxx = zeros(D*N);
invKxx = zeros(D*N);
Kww = zeros(D*M);
invKww = zeros(D*M);
for d=1:D
  % GP prior for X
  indDx = row2ind(D,N,d);
  Kxx(indDx,indDx) = kX{d}(inX,inX,thetaX{d});
  % GP prior for W
  indDw = col2ind(M,D,d);
  Kww(indDw,indDw) = kW{d}(inW,inW,thetaW{d});
end

% DEBUGGING: this helps to keep the covariance matrices numerically
% positive definite
reg = 1e-8
Kxx = Kxx + reg * eye(size(Kxx));
Kww = Kww + reg * eye(size(Kww));

cholKxx = safechol(Kxx);
cholKww = safechol(Kww);
logdetKxx = 2 * sum(log(diag(cholKxx)));
logdetKww = 2 * sum(log(diag(cholKww)));


%%
%% Initialize posterior parameters

X = zeros(D,N);
CovX = zeros(D*N);

W = randn(M,D);
CovW = zeros(M*D);

tau = 1;
logtau = 0;

mu = zeros(M,1);
v_mu = zeros(M,1);

logP = -inf;

%%
%% VB learning

for ind=1:opts.maxiter
  
  oldW = W;
  
  startitercpu = cputime;

  %% Update X
  sumKWWK = 0;
  sumKWY = 0;
  for n=1:N
    sumWW = 0;
    % Calculate <W'W> for a particular n
    for m=vecM(Obs(:,n))
      indMw = row2ind(M,D,m);
      sumWW = sumWW + W(m,:)'*W(m,:) + CovW(indMw,indMw);
    end
    indNx = col2ind(D,N,n);
    sumKWWK = sumKWWK + Kxx(:,indNx) * sumWW * Kxx(indNx,:);
    sumKWY = sumKWY + Kxx(:,indNx) * W(Obs(:,n),:)'*Y(Obs(:,n),n);
  end
%  [sqrtCovX, CovX, X(:)] = solvegpmvn(Kxx, tau*sumKWWK, tau * sumKWY);
  [sqrtCovX, CovX, X(:)] = solvegpmvn(Kxx, Kxx + tau*sumKWWK, tau * sumKWY);
  
  %% Update W
  sumKXXK = 0;
  sumKXY = 0;
  for m=1:M
    sumXX = 0;
    % Calculate <XX'> for a particular m
    for n=vecN(Obs(m,:))
      indNx = col2ind(D,N,n);
      sumXX = sumXX + X(:,n)*X(:,n)' + CovX(indNx,indNx);
    end
    indMw = row2ind(M,D,m);
    sumKXXK = sumKXXK + Kww(:,indMw) * sumXX * Kww(indMw,:);
    sumKXY = sumKXY + Kww(:,indMw) * X(:,Obs(m,:))*Y(m,Obs(m,:))';
  end
%  [sqrtCovW, CovW, W(:)] = solvegpmvn(Kww, tau*sumKXXK, tau * sumKXY);
  [sqrtCovW, CovW, W(:)] = solvegpmvn(Kww, Kww + tau*sumKXXK, tau * sumKXY);
  
  % CPU time used for calculations
  itercpu = cputime - startitercpu;
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Calculate lower bound of the log-likelihood

  startcostcpu = cputime;
  
  % Cost from Y
  cost_Y = 0;
  for n=1:N
    for m=vecM(Obs(:,n))
      mean_muY = W(m,:)*X(:,n) + mu(m);
      indNx = col2ind(D,N,n);
      indMw = row2ind(M,D,m);
      var_muY = W(m,:)*CovX(indNx,indNx)*W(m,:)' + ...
                X(:,n)'*CovW(indMw,indMw)*X(:,n) + ...
                traceprod(CovX(indNx,indNx), CovW(indMw,indMw)) + ...
                v_mu(m);
      tauY = tau;
      logtauY = logtau;
      cost_Y = cost_Y + mvnvbcost(Y(m,n),[], mean_muY,var_muY, tauY,logtauY);
    end
  end
    

  % Cost from W
  cost_W = mvnvbcost(W(:),CovW, 0,0, cholKww,logdetchol(cholKww));

  % Cost from X
  cost_X = mvnvbcost(X(:),CovX, 0,0, cholKxx,logdetchol(cholKxx));

  oldlogP = logP;

  % Lower bound of loglikelihood
  logP = cost_Y + cost_W + cost_X;% + cost_mu + cost_tau + cost_w;
  costcpu = cputime - startcostcpu;
  
  % Change in subspace
  angle = subspace(oldW, W);
  
  %% Print progress
  fprintf('%d step: loglike=%e, angle=%.2e (itercpu=%.2fs, nllcpu=%.2fs)\n', ...
          ind, logP, angle, itercpu, costcpu);
  
  if logP < oldlogP
    error('Cost increased!!');
  end
  
end

if nargout == 1
  Q.W = W;
  Q.CovW = CovW;
  Q.X = X;
  Q.CovX = CovX;
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

