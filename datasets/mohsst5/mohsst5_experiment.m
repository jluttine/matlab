function [Q_gppca, Q_vbpca] = mohsst5_experiment(data, D, maxiter, densityW, ...
                                                densityX, updatepseudow, ...
                                                updatepseudox)

% $$$ randn('state', 1)
% $$$ rand('state', 1)

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

datestring = datestr(now, 'yyyymmdd');

% Weight observations!
weights = sqrt(cosd(data.coordinates(2,:)));
Yobsw = bsxfun(@times, weights(:), Yobs);
Ytestw = bsxfun(@times, weights(:), Ytest);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% $$$ % Covariance functions
% $$$ % W distances in kilometers
% $$$ switch 1
% $$$  case 1
% $$$   gpcovW = {@gpcovScale, @(logtheta,x1,x2) gpcov(logtheta,x1,x2,@sqdistEarth)};
% $$$   logthetaW = [log(2); log(1e3)];
% $$$   Wcov = 'se';
% $$$ end
% $$$ % X distances in days
% $$$ switch 2
% $$$   case 2
% $$$    gpcovX = @gpcovPP;
% $$$    logthetaX = log(80);
% $$$    Xcov = 'cs';
% $$$ end
% $$$ 
% $$$ for d=1:D
% $$$   covfunW{d} = gpcovW;
% $$$   covfunX{d} = gpcovX;
% $$$   initthetaW{d} = logthetaW;
% $$$   initthetaX{d} = logthetaX;
% $$$ end
% $$$ 
% $$$ filename = sprintf(['/share/bayes/data/jaakko/mohsst5/mohsst5_testsize=%d_D=' ...
% $$$                     '%d_Wcov=%s_Wden=%d_Xcov=%s_Xden=%d_iter=%d_date=%s'], ...
% $$$                    testsize, D, Wcov, ceil(100*densityW), Xcov, ...
% $$$                    ceil(100*densityX), maxiter, datestring)
% $$$ 
% $$$ disp('-- Run GP PCA --')
% $$$ Q_gppca = vbgppcamv(Yobs, D,...
% $$$                     data.coordinates, data.time,...
% $$$                     covfunW, initthetaW,...
% $$$                     covfunX,initthetaX, ...
% $$$                     'maxiter',maxiter, ...
% $$$                     'pseudodensityx', densityX, ...
% $$$                     'pseudodensityw', densityW, ...
% $$$                     'loglikelihood', true, ...
% $$$                     'updatehyper', min(10,maxiter), ...
% $$$                     'updatepseudox', updatepseudox, ...
% $$$                     'updatepseudow', updatepseudow, ...
% $$$                     'maxsearchx', 3, ...
% $$$                     'maxsearchw', 3, ...
% $$$                     'checkgradx', false, ...
% $$$                     'autosavetime', 3600, ...
% $$$                     'autosavefile', filename);
% $$$ 
% $$$ 
% $$$ Yh_gppca = Q_gppca.W * Q_gppca.X;
% $$$ trainrmse_gppca = rmse(Yobs, Yh_gppca)
% $$$ testrmse_gppca = rmse(Ytest, Yh_gppca)
% $$$ % Huuuuuge matrices... Don't save them..
% $$$ Q_gppca.CovXp = [];
% $$$ Q_gppca.CovWp = [];
% $$$ 
% $$$ Q_gppca.filename = filename;
% $$$ 
% $$$ save(filename, '-struct', 'Q_gppca');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('-- Run PCA --')

% Alexander's vbpca
disp('version by Alex')
Q_vbpca = pca_full(Yobsw,D,'maxiters',maxiter,'rotate2pca',true, 'algorithm','vb');
Yh_vbpca = Q_vbpca.A * Q_vbpca.S + repmat(Q_vbpca.Mu,1,N);
vb_filename = sprintf(['/share/bayes/data/jaakko/mohsst5/' ...
                    'mohsst5_weighted_VBPCA_testsize=%d_D=%d_iter=%d_date=%s'], ...
                      testsize, D, maxiter, datestring);

% My vbpca
%disp('version by Jaakko')
%Q_vbpca = vbpcamv(Yobs,D,'maxiters',maxiter);
%Yh_vbpca = Q_vbpca.W * Q_vbpca.X + repmat(Q_vbpca.mu,1,N);

trainrmse_vbpca = rmse(Yobsw, Yh_vbpca)
testrmse_vbpca = rmse(Ytestw, Yh_vbpca)

% INVERSE WEIGHTS!!! Where is the inverese weighting for Av?????!!!!!
Q_vbpca.A = bsxfun(@rdivide, Q_vbpca.A, weights(:));
%Q_vbpca.varW = bsxfun(@rdivide, Q_vbpca.varW, weights(:).^2);
% $$$ for d=1:D
% $$$   lat = Q_gppca.pseudoW{d}(2,:);
% $$$   Q_gppca.Wp{d} = Q_gppca.Wp{d} ./ sqrt(abs(cosd(lat(:))));
% $$$ end

save(vb_filename, '-struct', 'Q_vbpca');
fprintf('Saved results in %s\n', vb_filename);