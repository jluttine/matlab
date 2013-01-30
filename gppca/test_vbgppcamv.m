function test_vbgppcamv(vbpca, ppca)

rand('state', 6);
randn('state', 6);

if nargin < 1
  vbpca = false;
end
if nargin < 2
  ppca = false;
end

inW = 1:100;
inX = 1:200;

M = length(inW);
N = length(inX);
D = 2;
X = zeros(D,N);
W = zeros(M,D);

% Covariance matrices
logthetaW = log(5);
logthetaX = log(1);

% Covariance functions
gpcovW = {@gpcovScale, @gpcov}; logthetaW = [3; logthetaW(:)];
%gpcovW = {@gpcovConstScale, 10, @gpcov};
%gpcovW = {@gpcov};
%gpcovX = {@gpcovScale, @gpcov}; logthetaX = [3; logthetaX(:)]
gpcovX = @gpcov;


% Generate latent variables
for d=1:D
% $$$   W(:,d) = mvnrnd(zeros(1,M),CovW)';
% $$$   X(d,:) = mvnrnd(zeros(1,N),CovX);
  lsW = 10*ceil(exp(logthetaW(end)));
  lsX = 30*ceil(exp(logthetaX(end)));
  W(:,d) = filter(hamming(lsW), 1, randn(M,1));
  X(d,:) = filter(hamming(lsX), 1, randn(N,1));
% $$$   W(:,d) = filter(ones(lsW,1)/lsW, 1, randn(M,1));
% $$$   X(d,:) = filter(ones(lsX,1)/lsX, 1, randn(N,1));
  covfunW{d} = gpcovW;
  covfunX{d} = gpcovX;
  initthetaW{d} = logthetaW-1;
  initthetaX{d} = logthetaX+log(2);
end

W = W * diag( sqrt(M./rowsum(W.^2)) );
X = diag( sqrt(N./colsum(X.^2)) ) * X;

% $$$ tsplot(X)
% $$$ tsplot(W')
% $$$ return

% Generate data
Y = W*X;

% Noise level
s2 = 0.1;

% Make noisy observations
Yn = Y + sqrt(s2)*randn(M,N);  % noise
Ynm = Yn; Ynm(rand(M,N)<0.5) = nan; % missing values
Ytest = Yn; Ytest(~isnan(Ynm)) = nan;

% Learn GP PCA
maxiter = 1e1;
Q = vbgppcamv(Ynm,D,inW,inX,covfunW,initthetaW,covfunX,initthetaX, ...
              'maxiter',maxiter, 'pseudodensityx', 0.9, 'pseudodensityw', ...
              0.4, 'loglikelihood', true, 'updatehyper', 5, 'maxsearchx', ...
              3, 'maxsearchw', 3, 'updatepseudow', true, 'updatepseudox', ...
              true, 'checkgradx', false);

est_noise = 1/Q.tau
real_noise = s2

% Learn other models
if ppca
  Qppca = pca_full(Ynm,D,'maxiters',maxiter,'rotate2pca',true, ...
                   'algorithm','ppca');
  Yh_ppca = Qppca.A * Qppca.S + repmat(Qppca.Mu,1,N);
  testrmse_ppca = rmse(Ytest, Yh_ppca)
  noiselessrmse_ppca = rmse(Y, Yh_ppca)
end
if vbpca
  Qvbpca = vbpcamv(Ynm,M-1,'maxiters',maxiter);
  Yh_vbpca = Qvbpca.W * Qvbpca.X + repmat(Qvbpca.mu,1,N);
  testrmse_vbpca = rmse(Ytest, Yh_vbpca)
  noiselessrmse_vbpca = rmse(Y, Yh_vbpca)
end



Yh_gppca = Q.W * Q.X;
testrmse_gppca = rmse(Ytest, Yh_gppca)
noiselessrmse_gppca = rmse(Y, Yh_gppca)
figure
subplot(5,1,1);
plot(Y');
title('Original noiseless data')
subplot(5,1,2);
plot(Ynm');
title('Observed data')
subplot(5,1,3);
plot(Yh_gppca');
title('Reconstruction of GP VB PCA')
if vbpca
  subplot(5,1,4);
  plot(Yh_vbpca');
  title('Reconstruction of VB PCA')
end
if ppca
  subplot(5,1,5);
  plot(Yh_ppca');
  title('Reconstruction of PPCA')
end

vX = Q.varX;%diag(Q.CovX);
eX = 2*sqrt(vX);%sqrt( vX(reshape(1:(N*D), D, N)) );
vW = Q.varW;%diag(Q.CovW);
eW = 2*sqrt(vW);%sqrt( vW(reshape(1:(M*D), M, D)) );
tsgpplot(inX, Q.X', eX', 'pseudoinputs', {Q.pseudoX});
tsgpplot(inW, Q.W, eW, 'pseudoinputs', {Q.pseudoW});

thetaX1 = exp(Q.logthetaX{1})
thetaX2 = exp(Q.logthetaX{2})

% $$$ exp(Q.logthetaW{1})
% $$$ exp(Q.logthetaW{2})
% $$$ 
% $$$ exp(Q.logthetaX{1})
% $$$ exp(Q.logthetaX{2})

% $$$ % DEBUG: (test the effect of CovXp)
% $$$ X = zeros(D,N);
% $$$ varX = zeros(D,N);
% $$$ for d=1:D
% $$$   [X(d,:), varX(d,:)] = gppred(Q.pseudoX{d}, Q.Xp{d}, 0, inX, covfunX{d}, ...
% $$$                                Q.logthetaX{d});
% $$$ end    
% $$$ varX;
% $$$ tsgpplot(inX, X', 2*sqrt(varX'), 'pseudoinputs', {Q.pseudoX});
% $$$ exp(Q.logthetaX{1})

% $$$ figure
% $$$ pcolor(Q.CovXp{1});
% $$$ figure
% $$$ pcolor(Q.CovWp{1});


return














%inv(funcK_noiseless(D,[10;10;1;10;10]));return
[p, loglike] = gplearn(ymv, @(p) funcK_noisy(D,p), [10;10;1;300;300;0.1])
%[p, loglike] = gplearn(ymv, @(p) funcK_noisy(D,p), [300;300;10;1;1])
%[p, loglike] = gplearn(ymv, @(p) funcK_noisy(D,p), [300;300;10;5;5;10;0.1])
%[p, loglike] = gplearn(ymv, @(p) funcK_noisy(D,p), [1000;100;0.1;3;100;0.059;1;0.1])
%[p, loglike] = gplearn(ymv, @(p) funcK_noisy(D,p), [100;100;3;100;0.059;1;0.1])

%p(1:4) = [2000;30;0.1;10]
%p(3) = 0.1;
[mu,Cov] = gppred(th,tmv,ymv, @(x1,x2) funcK_noiseless(gpdist(x1,x2),p(1:(end-1))), p(end));
yh = mnorm_rnd(mu, Cov, 10);

figure
clf
plot(th, yh, 'r')
figure
clf
gpplot(th,mu,Cov);
hold on
plot(tmv,ymv,'k+')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function gpmapplot(coord, y)
x = coord';
% Use haversine distance measure
mydist = @(x1,x2) gpdist(x1,x2,@dist_coord);
D = mydist(x,x);
funcK = @(p) gpK(@() gpK_ratquad(D,p(1),p(2),p(3)), ...
                 @() gpK_noise(length(D),p(4)));
[p, loglike] = gplearn(y, funcK, [10;5;1;1e-5])
funcK_noiseless = @(D) gpK(@() gpK_ratquad(D, p(1),p(2),p(3)));
[LONI, LATI] = get_grid(40, 40);
xh = [LONI(:)';LATI(:)'];
[mu,Cov] = gppred(xh,x,y, @(x1,x2) funcK_noiseless(mydist(x1,x2)), p(4));
yh = mu;
ZI = reshape(yh, size(LONI));
plot_map
hold on
m_pcolor(LONI, LATI, ZI);

% Set nice colormap
colormap(climcolmap);
shading flat;
lim = max( -min(yh), max(yh) );
set(gca, 'clim', [-lim lim]);

return

% $$$ %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% $$$ function D = gpdist(x1, x2, funcDist)
% $$$ if isvector(x1)
% $$$   x1 = x1(:)';
% $$$   x2 = x2(:)';
% $$$ end
% $$$ if nargin < 3
% $$$   funcDist = @(z1,z2) abs(z1-z2);
% $$$ end
% $$$ %[X2,X1] = meshgrid(x2,x1);
% $$$ n1 = size(x1,2);
% $$$ n2 = size(x2,2);
% $$$ D = zeros(n1,n2);
% $$$ for i=1:n1
% $$$   for j=1:n2
% $$$     D(i,j) = funcDist(x1(:,i),x2(:,j));
% $$$   end
% $$$ end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function d = dist_coord(coord1, coord2, varargin)
% d = dist(X1, X2)
% returns geographical distance in kilometers
% input vectors must be the same size and shape!!
% Approximates earth as a sphere.
% Distance is calculated for 

if nargin > 2
  nargin
end

% Quadratic mean radius from Wikipedia
R_avg = 6372.795477598;

% Convert to radians
q = pi / 180;
lon1 = coord1(1,:) * q;
lat1 = coord1(2,:) * q;
lon2 = coord2(1,:) * q;
lat2 = coord2(2,:) * q;

% Distance calculation (haversine law)
dlat = lat2 - lat1;
dlon = lon2 - lon1;
a = sin(dlat/2).^2 + cos(lat1).*cos(lat2).*(sin(dlon/2).^2);
c = 2 * atan2(sqrt(a), sqrt(1-a));
d = R_avg * c;

return

% $$$ %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% $$$ function [p,loglike] = gplearn(y, funcK, init_p)
% $$$ opts = optimset('GradObj', 'on');
% $$$ [p, negloglike] = fminunc(@(p) cost(y, funcK, p), init_p, opts);
% $$$ loglike = -negloglike;
% $$$ 
% $$$ %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% $$$ function [mu,Cov] = gppred(xh, x, y, funcK, noise)
% $$$ if nargin == 1
% $$$   x = [];
% $$$   y = [];
% $$$ end
% $$$ 
% $$$ Kxhx = funcK(xh,x);
% $$$ invKxx = inv(funcK(x,x)+noise^2*eye(length(x)));
% $$$ Kxhxh = funcK(xh,xh);
% $$$ if isempty(Kxhx)
% $$$   Kxhx = 0;
% $$$ end
% $$$ if isempty(invKxx)
% $$$   invKxx = 0;
% $$$ end
% $$$ mu = Kxhx*invKxx*y;
% $$$ if isempty(mu)
% $$$   mu = zeros(size(xh));
% $$$ end
% $$$ Cov = Kxhxh - Kxhx*invKxx*Kxhx';

% $$$ % $$$ %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% $$$ % $$$ function gpplot(x, mu, Cov)
% $$$ % $$$ e = sqrt(diag(Cov));
% $$$ % $$$ X = [x(1:end); x(end:-1:1)];
% $$$ % $$$ Y = [mu(1:end)+e; mu(end:-1:1)-e(end:-1:1)];
% $$$ % $$$ C = [.65,.65,.65];
% $$$ % $$$ fill(X,Y,C,'EdgeColor',C);
% $$$ % $$$ hold on
% $$$ % $$$ plot(x,mu,'k');
% $$$ 
% $$$ %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% $$$ function [K, dK] = gpK(varargin)
% $$$ K = 0;
% $$$ dK = [];
% $$$ for i=1:nargin
% $$$   if nargout == 1
% $$$     Knew = varargin{i}();
% $$$   else
% $$$     [Knew,dKnew] = varargin{i}();
% $$$     n = size(dKnew,3);
% $$$     if isempty(dK)
% $$$       dK = dKnew;
% $$$     else
% $$$       dK(:,:,end+(1:n)) = dKnew;
% $$$     end
% $$$   end
% $$$   K = K + Knew;
% $$$ end
% $$$ 
% $$$ %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% $$$ function [K, dK] = gpK_noise(n, p1);
% $$$ K = p1^2 * eye(n);
% $$$ if nargout >= 2
% $$$   dK = K * 2/p1;
% $$$ end
% $$$ 
% $$$ % $$$ %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% $$$ % $$$ function [K, dK] = gpK_sqexp(D, p1, p2)
% $$$ % $$$ K = p1^2*exp(-0.5*D.^2/(p2^2));
% $$$ % $$$ 
% $$$ % $$$ if nargout >= 2
% $$$ % $$$   dK = zeros([size(D), 2]);
% $$$ % $$$   dK(:,:,1) = K .* 2 / p1;
% $$$ % $$$   dK(:,:,2) = K .* (-0.5*D.^2) .* (-2*p2^(-3));
% $$$ % $$$ end
% $$$ 
% $$$ %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% $$$ function [K, dK] = gpK_decper(D, p1, p2, p3, p4)
% $$$ K = p1^2*exp(-0.5*(D.^2)/(p2^2)-2*sin(pi*D*p3).^2/(p4^2));
% $$$ 
% $$$ if nargout >= 2
% $$$   dK = zeros([size(D),4]);
% $$$   dK(:,:,1) = K .* 2 / p1;
% $$$   dK(:,:,2) = K .* (-0.5*D.^2) .* (-2*p2^(-3));
% $$$   dK(:,:,3) = K .* (-2*sin(pi*D*p3)*2) .* cos(pi*D*p3) .* (pi*D);
% $$$   dK(:,:,4) = K .* (-2*sin(pi*D*p3).^2) * (-2)*p4^(-3);
% $$$ end
% $$$ 
% $$$ % $$$ %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% $$$ % $$$ function [K,dK] = gpK_ratquad(D, p1, p2, p3)
% $$$ % $$$ %%% p3 = p3^2;
% $$$ % $$$ f = 1 + D.^2/(2*p3^2*p2^2);
% $$$ % $$$ K = p1^2*f.^(-p3^2);
% $$$ % $$$ 
% $$$ % $$$ if nargout >= 2
% $$$ % $$$   dK = zeros([size(D),3]);
% $$$ % $$$   dK(:,:,1) = K .* 2 / p1;
% $$$ % $$$   dK(:,:,2) = K .* (-p3^2).*f.^(-1) .* (-2).*D.^2/(2*p3^2)*p2^(-3);
% $$$ % $$$   dK(:,:,3) = K .* (-2*p3*log(f) + (-p3^2)./f.*D.^2/(2*p2^2) * (-2) * p3^(-3));
% $$$ % $$$ end
% $$$ 
% $$$ %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% $$$ function [f, df] = cost(y,funcK,p)
% $$$ [K,dK] = funcK(p);
% $$$ [f,df] = loglikelihood(y,K,dK);
% $$$ f = -f;
% $$$ df = -df;
% $$$ 
% $$$ %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% $$$ function [loglike, dloglike] = loglikelihood(y,K,dK)
% $$$ invK = inv(K);
% $$$ logdetK = log(det(K));
% $$$ if logdetK < -1e100;
% $$$   logdetK = -1e100;
% $$$ end
% $$$ n = length(y);
% $$$ 
% $$$ loglike = -0.5*y'*invK*y - 0.5*logdetK - 0.5*n*log(2*pi);
% $$$ 
% $$$ if nargout >= 2
% $$$   m = size(dK,3);
% $$$   dloglike = zeros(m,1);
% $$$   a = invK * y;
% $$$   W = a*a' - invK;
% $$$   for i = 1:m
% $$$     dloglike(i) = 0.5*sum(sum(W.*dK(:,:,i)));
% $$$   end
% $$$ end
% $$$ 
% $$$ 
% $$$ % $$$ %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% $$$ % $$$ function y = mnorm_rnd(mu, Cov, n)
% $$$ % $$$ m = max(length(mu), length(Cov));
% $$$ % $$$ 
% $$$ % $$$ opts.issym = true;
% $$$ % $$$ opts.isreal = true;
% $$$ % $$$ [V,D] = svd(Cov);
% $$$ % $$$ D(D<0) = 0;
% $$$ % $$$ D = sqrt(D);
% $$$ % $$$ A = V * D;
% $$$ % $$$ 
% $$$ % $$$ y = repmat(mu,1,n) + A*randn(m,n);
