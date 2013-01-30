

% Test variational RPPCA by Jaakko

% TODO: This kind of test function for which the user can interactively
% give some parameters (N,M,D,ESTIMATED D,SEED,...)

% TODO: Similar interactive function for using VBRFA for real problems?
% It just asks the name of the data matrix and some other relevant
% parameters and then runs.

function test_vbrfa

n = 400;
m = 40;
d = 10;

randn('state', 2);
rand('state', 2);

mu = 5 * randn(m,1);
%W = randn(m,d) * diag(1*ones(d,1)); %diag(10:-1:(10-d+1));
%X = randn(d,n);

% Setup noiseless covariance matrix (singular values 2^2-1, 3^2-1,...)
lambda = zeros(m,1);
lambda(1:d) = 2^2 %((d:-1:1)' + 1) .^ 2 - 1;
V = orth(randn(m));
Cov = V*diag(lambda)*V';

% Generate normal noiseless data (sing.values 1^2, ..., 1^2, 2^2, 3^2, 4^2, ...)
%Y = bsxfun(@plus, W*X, mu);
Y = mvnrnd(mu(:)', Cov, n)';
varY = var(Y,1,2)

% Add unit variance noise
Yn = Y + randn(m,n);

% Generate some outliers
Yno = Yn;
p = (rand(m,n) < 0.01);
signs = ones(sum(p(:)),1);
signs(rand(sum(p(:)),1)<0.5) = -1;
Yno(p) = Yno(p) + 10*signs;
% $$$ Yno(p) = signs .* ( 50 + 20*randn(sum(p(:)),1) );

% Generate some missing values
Ynom = Yno;
pmv = (rand(m,n) < 0.1);
Ynom(pmv) = NaN;

% Use STANDARD PCA
[W_pca, X_pca, mu_pca, s2_pca] = pca_full(Yn, d);
Y_pca = W_pca*X_pca + repmat(mu_pca, 1, n);

plot_outliermode = false;

% Use ROBUST PROBABILISTIC PCA

%
% THIS INITIALIZATION SEEMS TO BE IMPORTANT?
%
init.tau = 10;
init.nu = 0.1;

% TODO: ROTATION DOES NOT SEEM TO BE VERY GOOD IN PRACTICE..? IS THERE
% SOME PROBLEMS IN IT? OK, NOW IT SEEMS TO BE WORKING..
%
% ROTATING TOO EARLY MAY LEAD TO TOO EFFICIENT PRUNING OUT OF
% COMPONENTS... WAITING EVEN FOR 1-2 ITERATIONS AT THE BEGINNING MAY BE
% SUFFICIENT. HOWEVER, IF MU IS INITIALIZED VERY BADLY, IT WOULD BE
% NECESSARY (MAYBE?) TO DO THE BIAS MOVING..?
%
% OR: INIT WITH LARGE NOISE VARIANCE AND ROTATING IMMEDIATELY LEADS TO
% PRUNING OUT COMPONENTS TOO MUCH..?
Q = vbrfa(Ynom, m-1, 'maxiter', 50, 'update_nu', 5, 'init', init, ...
          'rotate', 5);
% Q = vbrfa(Ynom, m-1, 'maxiters', 50, 'prior', prior, 'init', init);
Wh = Q.W;
Xh = Q.X;
Sh = covarray_to_covcell(Q.CovX);
muh = Q.Mu
nuh = Q.nu
alphah = Q.alpha;
sorted_alphah = sort(Q.alpha)
%s2h = 1./results.bTau;
%u = results.U;
Yh = Wh*Xh + repmat(muh,1,n);

% $$$ mean_X = mean(Xh)
% $$$ WW = Wh'*Wh

Wh;
muh;
%s2h;
nuh;
Sh;

[tmp,I] = sort(alphah,'ascend');
err_W = 180 * subspace(Wh(:,I(1:d)), W_pca) / pi
err_mu = sqrt(mean((muh-mu_pca).^2));

% $$$     % Plot using PCA
% $$$     if plot_outliermode
% $$$       V = [eye(2); zeros(m-2, 2)];
% $$$     else
% $$$       [V,tmp,D] = princomp(Yn');
% $$$       [D,I] = sort(D, 'descend');
% $$$       V = V(:,I(1:2));
% $$$     end
% $$$ 
% $$$     % Show in the real space
% $$$     VYh = V'*Yh;
% $$$     for i=1:length(Sh)
% $$$       VYvh{i} = V' * Wh * Sh{i} * Wh' * V;
% $$$     end
% $$$     if plot_outliermode
% $$$       VYno = V'*Yno;
% $$$     %subspace2d(VYh, VYvh, VYno);
% $$$     %xl = xlim; yl = ylim;
% $$$     %VYno(1,pmv(1,:)) = xl(end); % draw missing dimensions to borders
% $$$     %VYno(2,pmv(2,:)) = yl(end);
% $$$     %close gcf;
% $$$     else
% $$$       VYno = V'*Y;
% $$$     end
% $$$     subspace2d(VYh, VYvh, VYno);
% $$$     set(gca, 'DataAspectRatioMode', 'manual');
% $$$     set(gca, 'DataAspectRatio', [1 1 1]);
% $$$     % Mark outliers
% $$$     hold on
% $$$     indeces = sum(p,1)>0;
% $$$     plot(VYno(1,indeces), VYno(2, indeces), 'go', 'MarkerSize', 8);
% $$$ 
% $$$     % Mark observations with missing values
% $$$     indeces = sum(pmv,1)>0;
% $$$     scatter(VYno(1,indeces), VYno(2, indeces), 'yo');
% $$$ 
% $$$     title(sprintf('My RVBPCA n=%d m=%d d=%d outliers=%.3f mv=%.3f', n,m,d, ...
% $$$                   sum(p(:))/(n*m), sum(pmv(:))/(n*m)));

err_Ymy = sqrt(mean((Y(:)-Yh(:)).^2))
err_Ypca = sqrt(mean((Y(:)-Y_pca(:)).^2))

tsplot(Y,'k');
addtsplot(Ynom,'g');
addtsplot(Yh,'r');
%addtsplot(Y_pca,'g');

%tsplot(Y-Yh, 'k')

return

