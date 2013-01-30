function Q_gppca = mohsst5_experiment_vbinit(data, D, maxiter, densityW, ...
                                             densityX, updatepseudow, ...
                                             updatepseudox)
if nargin < 1
  maxiter = 10;
end
if nargin < 6
  updatepseudow = true;
end
if nargin < 7
  updatepseudox = true;
end

[M,N] = size(data.observations);

% Divide data into test and train sets
testsize = 20; % testset size in percents
Itest = load(sprintf('/share/bayes/data/jaakko/mohsst5/ind%dtest.mat',testsize));
Itrain = load(sprintf('/share/bayes/data/jaakko/mohsst5/ind%dtrain.mat',testsize));
Ytest = data.observations;
Ytest(Itrain.Itrain) = nan;
Yobs = data.observations;
Yobs(Itest.Itest) = nan;

% Run VB PCA
disp('-- Run VB PCA --')
Q_vbpca = pca_full(Yobs,D,'maxiters',min(10,maxiter),'rotate2pca',true, ...
                   'algorithm','vb');


% $$$ % Load VB PCA results with 80 components
% $$$ Q_vbpca = load('/home/alexilin/matlab/kaplan/vbresults80.mat');

disp('Initialize with VB PCA results');
Qinit.W = Q_vbpca.A(:,1:D);
Qinit.X = Q_vbpca.S(1:D,:);
Qinit.tau = 1/Q_vbpca.V;
Qinit.varW = zeros(size(Qinit.W));
Qinit.varX = zeros(size(Qinit.X));
for m=1:M
  v = diag(Q_vbpca.Av{m});
  Qinit.varW(m,:) = v(1:D);
end
for n=1:N
  v = diag(Q_vbpca.Sv{n});
  Qinit.varX(:,n) = v(1:D);
end
% Put the mean as the last component.. :|
Qinit.X(end,:) = 1;
Qinit.varX(end,:) = 0;
Qinit.W(:,end) = Q_vbpca.Mu;
Qinit.varW(:,end) = Q_vbpca.Muv;

% Covariance functions
% W distances in kilometers
switch 1
 case 1
  gpcovW = {@gpcovScale, @(logtheta,x1,x2) gpcov(logtheta,x1,x2,@sqdistEarth)};
  logthetaW = [log(2); log(1e3)];
  Wcov = 'se';
end
% X distances in days
switch 2
  case 2
   gpcovX = @gpcovPP;
   logthetaX = log(80);
   Xcov = 'cs';
end

for d=1:D
  covfunW{d} = gpcovW;
  covfunX{d} = gpcovX;
  initthetaW{d} = logthetaW;
  initthetaX{d} = logthetaX;
end

datestring = datestr(now, 'yyyymmdd');
filename = sprintf(['/share/bayes/data/jaakko/mohsst5/' ...
                    'mohsst5_vbinit_testsize=%d_D=%d_Wcov=%s_Wden=' ...
                    '%d_Xcov=%s_Xden=%d_iter=%d_date=%s'], testsize, D, ...
                   Wcov, ceil(100*densityW), Xcov, ceil(100*densityX), ...
                   maxiter, datestring)

disp('-- Run GP PCA --')
Q_gppca = vbgppcamv(Yobs, D,...
                    data.coordinates, data.time,...
                    covfunW, initthetaW,...
                    covfunX,initthetaX, ...
                    'init', Qinit, ...
                    'maxiter',maxiter, ...
                    'pseudodensityx', densityX, ...
                    'pseudodensityw', densityW, ...
                    'loglikelihood', true, ...
                    'updatehyper', [1 1 10:10:(maxiter-1)], ...
                    'updatepseudox', updatepseudox, ...
                    'updatepseudow', updatepseudow, ...
                    'maxsearchx', 3, ...
                    'maxsearchw', 3, ...
                    'checkgradx', false, ...
                    'autosavetime', 3600, ...
                    'autosavefile', filename);


Yh_gppca = Q_gppca.W * Q_gppca.X;
%Yh_vbpca = Q_vbpca.A * Q_vbpca.S + repmat(Q_vbpca.Mu,1,N);

trainrmse_gppca = rmse(Yobs, Yh_gppca)
%trainrmse_vbpca = rmse(Yobs, Yh_vbpca)
testrmse_gppca = rmse(Ytest, Yh_gppca)
%testrmse_vbpca = rmse(Ytest, Yh_vbpca)

% Huuuuuge matrices... Don't save them..
Q_gppca.CovXp = [];
Q_gppca.CovWp = [];

Q_gppca.filename = filename;

%save(filename, '-struct', 'Q_gppca');
save(filename, '-struct', 'Q_gppca');