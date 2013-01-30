function [Q_gppca, Q_vbpca] = mohsst5_experiment_periodic(data, D, maxiter, ...
                                                  densityW, densityX, ...
                                                  updatepseudow, updatepseudox, ...
                                                  vbinit, covfuncset)

seed = ceil(sum(clock))
randn('state', seed)
rand('state', seed)

if nargin < 1
  maxiter = 10;
end
if nargin < 6
  updatepseudow = true;
end
if nargin < 7
  updatepseudox = true;
end
if nargin < 8
  vbinit = false;
end
if nargin < 9
  covfuncset = 3;
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

init = [];
switch vbinit
 case 0
  vbinitstring = '';
 case 1
  Qvb = pca_full(Yobsw,D,'maxiters',10,'rotate2pca',true, 'algorithm','vb');
  init.W = zeros(M,D);
  init.X = randn(D,N);
  init.W(:,1:10) = Qvb.A(:,1:10);
  init.X(1:10,:) = Qvb.S(1:10,:);
  vbinitstring = '_vbinit';
 case 2
  Qvb = pca_full(Yobsw,D,'maxiters',10,'rotate2pca',true, 'algorithm','vb');

  init.W = zeros(M,D);
  init.X = randn(D,N);
  
  samplerate = 1/30; % one sample in month
  
  % 3 very slow components
  lowpass = 1 ./ (10*365);
  [Af,Sf,S] = dss(Qvb.A, Qvb.S, Qvb.Sv, 2 * lowpass/samplerate, 3);
  init.W(:,1:3) = Af;
  init.X(1:3,:) = S;

  % 5 slow components
  lowpass = 1 ./ [10*365, 2*365];
  [Af,Sf,S] = dss(Qvb.A, Qvb.S, Qvb.Sv, 2 * lowpass/samplerate, 5);
  init.W(:,4:8) = Af;
  init.X(4:8,:) = S;
  
  % 2 periodical components
  bandpass = 1 ./ [2*365, 0.5*365];
  [Af,Sf,S] = dss(Qvb.A, Qvb.S, Qvb.Sv, 2 * bandpass/samplerate, 2);
  init.W(:,9:10) = Af;
  init.X(9:10,:) = S;
  
% $$$   % DEBUGGING STUFF
% $$$   Q_gppca = Af;%init.W;
% $$$   Q_vbpca.Xf = Sf;
% $$$   Q_vbpca.X = S;
% $$$   return
% $$$   
  vbinitstring = '_vbdssinit';
  
end

% Covariance functions
% W distances in kilometers
for d=1:D
  switch 1 % select covariance functions
   case 1
    gpcovW = {@gpcovScale, @(logtheta,x1,x2) gpcov(logtheta,x1,x2,@ ...
                                                   sqdistEarth)};
    logthetaW = [log(2); log(1e3)];
    Wcov = 'se';
  end
  initthetaW{d} = logthetaW;
  covfunW{d} = gpcovW;
end

% X distances in days
for d=1:D
  switch covfuncset % select covariance functions
   case 1
    % 80 pediodical
    if d==1
      disp('Using 80 periodical');
    end
    gpcovX = {@gpcovProduct, @gpcov, @gpcovPeriodic};
    logthetaX = log([1e5, 0.5, 365*(D-d+3)]);
    Xcov = 'per';
   case 2
    % 80 piecewise polynomial
    if d==1
      disp('Using 80 piecewise polynomial');
    end
    gpcovX = @gpcovPP;
    logthetaX = log(30*(5+40*rand));
    Xcov = 'cs';
   case 3
    % 5 trend SE + 5 slow SE + 5 periodical + rest faster piecewise
    % polynomial
    if d==1
      disp(['Using 5 trends + 5 slow + 5 periodical + rest faster ' ...
            'piecewise polynomial']);
    end
    if d <= 5
      gpcovX = @gpcov;
      logthetaX = log(365*(10+10*rand)); % 10-20 years
      densityX(d) = 0.05;
    elseif d <= 10
      gpcovX = @gpcov;
      logthetaX = log(365*(2+8*rand)); % 2-10 years
      densityX(d) = 0.2;
    elseif d <= 15
      gpcovX = {@gpcovProduct, @gpcov, @gpcovPeriodic};
      logthetaX = log([200*365, 0.5, 365*(0.5+1.5*rand)]); % 0.5-2 years
      densityX(d) = 0.2;
    else
      gpcovX = @gpcovPP;
      logthetaX = log(30*(6+6*rand)); % 0.5-1 years
      densityX(d) = 1;
    end
    Xcov = 'se-per-cs';
   case 4
    % 3 slow piecewise polynomial + 5 middle + 2 periodical + rest faster
    % piecewise polynomial
    if d==1
      disp(['Using 3 trends + 5 slow + 2 periodical + rest faster ' ...
            'piecewise polynomial']);
    end
    if d <= 3
      gpcovX = @gpcov;
      logthetaX = log(365*(10+10*rand)); % 10-20 years
      densityX(d) = 0.05;
    elseif d <= 8
      gpcovX = @gpcov;
      logthetaX = log(365*(2+8*rand)); % 2-10 years
      densityX(d) = 0.2;
    elseif d <= 10
      gpcovX = {@gpcovProduct, @gpcov, @gpcovPeriodic};
      logthetaX = log([80*365, 0.5, 365*(0.5+1.5*rand)]); % 0.5-2 years
      densityX(d) = 0.2;
    else
      gpcovX = @gpcovPP;
      logthetaX = log(30*(6+6*rand)); % 0.5-1 years
      densityX(d) = 1;
    end
    Xcov = 'se-per-cs';
  end
  initthetaX{d} = logthetaX;
  covfunX{d} = gpcovX;
end


filename = sprintf(['/share/bayes/data/jaakko/mohsst5/' ...
                    'mohsst5_weighted%s_testsize=%d_D=%d_Wcov=%s_Wden=' ...
                    '%d_Xcov=%s_Xden=%d_iter=%d_date=%s'], vbinitstring, ...
                   testsize, D, Wcov, ceil(100*min(densityW)), Xcov, ...
                   ceil(100*min(densityX)), maxiter, datestring)

disp('-- Run GP PCA --')
Q_gppca = vbgppcamv(Yobsw, D,...
                    data.coordinates, data.time,...
                    covfunW, initthetaW,...
                    covfunX,initthetaX, ...
                    'maxiter',maxiter, ...
                    'pseudodensityx', densityX, ...
                    'pseudodensityw', densityW, ...
                    'loglikelihood', true, ...
                    'updatehyper', [1 6 21 51 101 201 501], ...
                    'updatepseudox', updatepseudox, ...
                    'updatepseudow', updatepseudow, ...
                    'maxsearchx', 3, ...
                    'maxsearchw', 3, ...
                    'checkgradx', false, ...
                    'checkgradw', false, ...
                    'init', init, ...
                    'autosavetime', 3600, ...
                    'autosavefile', filename);


% Error measures (using weights!)
Yh_gppca = Q_gppca.W * Q_gppca.X;
trainrmse_gppca = rmse(Yobsw, Yh_gppca)
testrmse_gppca = rmse(Ytestw, Yh_gppca)

% Huuuuuge matrices... Don't save them..
Q_gppca.CovXp = [];
Q_gppca.CovWp = [];

% INVERSE WEIGHTS!!!
Q_gppca.W = bsxfun(@rdivide, Q_gppca.W, weights(:));
Q_gppca.varW = bsxfun(@rdivide, Q_gppca.varW, weights(:).^2);
for d=1:D
  lat = Q_gppca.pseudoW{d}(2,:);
  Q_gppca.Wp{d} = Q_gppca.Wp{d} ./ sqrt(abs(cosd(lat(:))));
end

% Save
Q_gppca.filename = filename;
Q_gppca.seed = seed;
save(filename, '-struct', 'Q_gppca');
fprintf('Saved results in %s\n', filename);

% $$$ %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% $$$ disp('-- Run PCA --')
% $$$ 
% $$$ % Alexander's vbpca
% $$$ disp('version by Alex')
% $$$ Q_vbpca = pca_full(Yobs,D,'maxiters',maxiter,'rotate2pca',true, 
%'algorithm','vb');
% $$$ Yh_vbpca = Q_vbpca.A * Q_vbpca.S + repmat(Q_vbpca.Mu,1,N);
% $$$ vb_filename = sprintf(['/share/bayes/data/jaakko/mohsst5/' ...
% $$$                     
%'mohsst5_VBPCA_testsize=%d_D=%d_iter=%d_date=%s'], ...
% $$$                       testsize, D, maxiter, datestring);
% $$$ save(vb_filename, '-struct', 'Q_vbpca');
% $$$ 
% $$$ % My vbpca
% $$$ %disp('version by Jaakko')
% $$$ %Q_vbpca = vbpcamv(Yobs,D,'maxiters',maxiter);
% $$$ %Yh_vbpca = Q_vbpca.W * Q_vbpca.X + repmat(Q_vbpca.mu,1,N);
% $$$ 
% $$$ trainrmse_vbpca = rmse(Yobs, Yh_vbpca)
% $$$ testrmse_vbpca = rmse(Ytest, Yh_vbpca)


