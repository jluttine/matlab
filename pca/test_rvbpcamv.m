

% Test variational RPPCA by Jaakko

n = 400;
m = 5;
d = 3;

randn('state', 1);
rand('state', 1);

mu = 20 * randn(m,1);
B = randn(m,d);
W = orth(randn(m,d)) * diag(10:-1:(10-d+1));
X = orth(randn(d,n)')' * sqrt(n);
%X = X - repmat(mean(X,2),1,n);
%X = X ./ repmat(std(X',1)',1,n);
s2 = 0.1;

% Generate normal data
Y = W*X + repmat(mu,1,n);

% Add noise
Yn = Y + sqrt(s2)*randn(m,n);

% Generate some outliers
Yno = Yn;
p = (rand(m,n) < 0.01);
Yno(p) = -200 * rand(sum(p(:)),1);

% Generate some missing values (but don't lose outliers :))
Ynom = Yno;
pmv = (rand(m,n) < 0.1);% & (~p);
Ynom(pmv) = NaN;

% Use STANDARD PCA
[W_pca, X_pca, mu_pca, s2_pca] = pca_full(Yn, d);
Y_pca = W_pca*X_pca + repmat(mu_pca, 1, n);

plot_outliermode = false;

% Use ROBUST PROBABILISTIC PCA
    prior = [];
    Q = rvbpcamv(Ynom, m-1, 'maxiters', 100, 'prior', prior);
    Wh = Q.W;
    Xh = Q.X;
    Sh = Q.Sv;
    muh = Q.mu;
    nuh = Q.nu;
    %s2h = 1./results.bTau;
    %u = results.U;
    Yh = Wh*Xh + repmat(muh,1,n);

    Wh;
    muh;
    %s2h;
    nuh;
    Sh;

    err_W = 180 * subspace(Wh(:,1:d), W_pca) / pi
    err_mu = sqrt(mean((muh-mu_pca).^2));

    % Plot using PCA
    if plot_outliermode
      V = [eye(2); zeros(m-2, 2)];
    else
      [V,tmp,D] = princomp(Yn');
      [D,I] = sort(D, 'descend');
      V = V(:,I(1:2));
    end

    % Show in the real space
    VYh = V'*Yh;
    for i=1:length(Sh)
      VYvh{i} = V' * Wh * Sh{i} * Wh' * V;
    end
    if plot_outliermode
      VYno = V'*Yno;
    %subspace2d(VYh, VYvh, VYno);
    %xl = xlim; yl = ylim;
    %VYno(1,pmv(1,:)) = xl(end); % draw missing dimensions to borders
    %VYno(2,pmv(2,:)) = yl(end);
    %close gcf;
    else
      VYno = V'*Y;
    end
    subspace2d(VYh, VYvh, VYno);
    set(gca, 'DataAspectRatioMode', 'manual');
    set(gca, 'DataAspectRatio', [1 1 1]);
    % Mark outliers
    hold on
    indeces = sum(p,1)>0;
    plot(VYno(1,indeces), VYno(2, indeces), 'go', 'MarkerSize', 8);

    % Mark observations with missing values
    indeces = sum(pmv,1)>0;
    scatter(VYno(1,indeces), VYno(2, indeces), 'yo');

    title(sprintf('My RVBPCA n=%d m=%d d=%d outliers=%.3f mv=%.3f', n,m,d, ...
                  sum(p(:))/(n*m), sum(pmv(:))/(n*m)));
    
    err_Ymy = sqrt(mean((Y(:)-Yh(:)).^2))
    err_Ypca = sqrt(mean((Y(:)-Y_pca(:)).^2))

return

