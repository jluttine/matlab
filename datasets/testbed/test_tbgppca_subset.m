function test_tbgppca_subset

data = tbload_temperature;
data = tbpreprocess(data);
data = tbsubset(data, 1:1:79, 1:10:89000);

data

% $$$ figure
% $$$ plot(data.observations');

inW = data.coordinates(:,1:2)';
inX = data.time';

[M,N] = size(data.observations);
D = 2;

% Covariance matrices
logthetaW = {[log(1);log(5e1)], [log(1);log(1e1)]}
logthetaX = {log(2000), log(100)};

% Covariance functions
gpcovW = {@gpcovScale, @(logtheta,x1,x2) gpcov(logtheta,x1,x2,@sqdistEarth)};
% $$$ pseudoW = cell(D,1);%zeros([dimW,Mp,D]);
% $$$ Mp = ceil(1*M);
% $$$ for d=1:D
% $$$   permM = randperm(M);
% $$$   pseudoW{d} = inW(:,permM(1:Mp));% + randn(2,Mp);
% $$$ end

gpcovX = @gpcov;

% Generate latent variables
for d=1:D
  covfunW{d} = gpcovW;
  covfunX{d} = gpcovX;
end

Y = data.observations;

% Learn GP PCA
maxiter = 3e1;
Q = vbgppcamv(Y,D,inW,inX,covfunW,logthetaW,covfunX,logthetaX, ...
              'maxiter',maxiter, 'pseudodensityx', 0.01, 'pseudodensityw', ...
              0.4, 'loglikelihood', true, 'updatehyper', 5, 'maxsearchx', ...
              3, 'maxsearchw', 50);

% $$$ % Learn other models
% $$$ if ppca
% $$$   Qppca = pca_full(Y,D,'maxiters',maxiter,'rotate2pca',true, ...
% $$$                    'algorithm','ppca');
% $$$   Yh_ppca = Qppca.A * Qppca.S + repmat(Qppca.Mu,1,N);
% $$$   testrmse_ppca = rmse(Ytest, Yh_ppca)
% $$$   noiselessrmse_ppca = rmse(Y, Yh_ppca)
% $$$ end
% $$$ if vbpca
% $$$   Qvbpca = vbpcamv(Ynm,M-1,'maxiters',maxiter);
% $$$   Yh_vbpca = Qvbpca.W * Qvbpca.X + repmat(Qvbpca.mu,1,N);
% $$$   testrmse_vbpca = rmse(Ytest, Yh_vbpca)
% $$$   noiselessrmse_vbpca = rmse(Y, Yh_vbpca)
% $$$ end



Yh_gppca = Q.W * Q.X;
rmse_gppca = rmse(Y, Yh_gppca)
%testrmse_gppca = rmse(Ytest, Yh_gppca)
%noiselessrmse_gppca = rmse(Y, Yh_gppca)

figure
subplot(5,1,1);
plot(Y');
title('Original data')
% $$$ subplot(5,1,2);
% $$$ plot(Ynm');
% $$$ title('Observed data')
subplot(5,1,3);
plot(Yh_gppca');
title('Reconstruction of GP VB PCA')
% $$$ if vbpca
% $$$   subplot(5,1,4);
% $$$   plot(Yh_vbpca');
% $$$   title('Reconstruction of VB PCA')
% $$$ end
% $$$ if ppca
% $$$   subplot(5,1,5);
% $$$   plot(Yh_ppca');
% $$$   title('Reconstruction of PPCA')
% $$$ end

vX = Q.varX;%diag(Q.CovX);
eX = 2*sqrt(vX);%sqrt( vX(reshape(1:(N*D), D, N)) );
vW = Q.varW;%diag(Q.CovW);
eW = 2*sqrt(vW);%sqrt( vW(reshape(1:(M*D), M, D)) );
tsgpplot(inX, Q.X', eX', 'pseudoinputs', {Q.pseudoX});

%tsgpplot(inW, Q.W, eW, 'pseudoinputs', {Q.pseudoW});
figure
for d=1:D
  subplot(D,1,d);
  mapproj('testbed');
  mapplot(inW(1,:),inW(2,:),'ro');
  hold on
  mapplot(Q.pseudoW{d}(1,:),Q.pseudoW{d}(2,:),'k+');
  mapcoast
end

figure
mapproj('testbed');
gpmapcolor(Q.pseudoW, Q.Wp, Q.CovWp, 22:0.1:27, 59:0.1:61, Q.logthetaW, covfunW);
