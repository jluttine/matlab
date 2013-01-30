function loglike = gppca_loglikelihood(Y, Qgp, Ytrain)
% loglike = gppca_loglikelihood(Y, Qgp)
%
% Returns posterior predictive loglikelihood.

% Get covariance matrices (needs the training data)
D = cols(Qgp.W);
if isempty(Qgp.CovWp) || isempty(Qgp.CovXp)
  disp('Must run GP PCA to get covariance matrices.');
  Q = vbgppcamv(Ytrain,D,Qgp.inW,Qgp.inX,Qgp.covfuncW,Qgp.logthetaW,Qgp.covfuncX, ...
                Qgp.logthetaX,'init',Qgp, 'maxiter',1, 'updatepseudox', ...
                false, 'updatepseudow', false, 'initpseudow', {Qgp.pseudoW}, ...
                'initpseudox', {Qgp.pseudoX}, 'loglikelihood', false, ...
                'updatehyper', false, 'reconstruct', false);
  Qgp.CovXp = Q.CovXp;
  Qgp.CovWp = Q.CovWp;  
end

% Remove empty rows/columns
Obs = ~isnan(Y);
rowrm = (colsum(Obs) == 0);
colrm = (rowsum(Obs) == 0);
if any(rowrm)
  Qgp.inW(:,rowrm) = [];
  Y(rowrm,:) = [];
  Obs(rowrm,:) = [];
end
if any(colrm)
  Qgp.inX(:,colrm) = [];
  Y(:,colrm) = [];
  Obs(:,colrm) = [];
end

% Sample X and W
N = 10; % number of samples
W = zeros([rows(Y), D, N]);
X = zeros([D, cols(Y), N]);
for d=1:cols(Qgp.W)
  W(:,d,:) = reshape(gppredrnd(Qgp.inW, Qgp.pseudoW{d}, Qgp.Wp{d}, ...
                               Qgp.CovWp{d}, Qgp.logthetaW{d}, Qgp.covfuncW{d}, ...
                               N), [rows(W), 1, N]);
  X(d,:,:) = reshape(gppredrnd(Qgp.inX, Qgp.pseudoX{d}, Qgp.Xp{d}, ...
                               Qgp.CovXp{d}, Qgp.logthetaX{d}, Qgp.covfuncX{d}, ...
                               N), [1, cols(X), N]);
  
  fprintf('Sampled dimension %d/%d\n', d, D);
end

% Sample tau
tau = gamrnd(Qgp.a_tau, 1/Qgp.b_tau, [N,1]);

% Estimate predictive loglikelihood by averaging over samples
loglikes = zeros(N,1);
loglike = 0;
[I,J] = find(Obs);
for k=1:length(I)
  for n=1:N
    mu = W(I(k),:,n)*X(:,J(k),n); % reconstruct
    loglikes(n) = sum(norm_lpdf(Y(I(k),J(k)), mu, sqrt(1/tau(n))));
  end
  loglike = loglike + log( mean(exp(loglikes)) );
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function X = gppredrnd(inputpred, input, postmu, postCov, logtheta, covfunc, N)

I = eye(length(postmu));

if ~iscell(covfunc), covfunc = {covfunc}; end

% Calculate the covariance matrices
Kpp = feval(covfunc{:}, logtheta, input, input);
Kpp = Kpp + 1e-6*I;
Lp = chol(Kpp, 'lower');

Kxp = feval(covfunc{:}, logtheta, inputpred, input);

Kxx = feval(covfunc{:}, logtheta, inputpred, inputpred);


R = solve_tril(Lp, Kxp');
S = solve_triu(Lp', R);

mu = R' * solve_tril(Lp, postmu);
% $$$ szKxx = size(Kxx)
% $$$ szS = size(S)
% $$$ szKpp = size(Kpp)
% $$$ szPostCov = size(postCov)
Cov = Kxx - S' * (Kpp - postCov) * S + 1e-6*eye(length(mu));

X = zeros(length(Kxx), N);
for i=1:N
  X(:,i) = mymvnrnd(mu, Cov);
end
