
function Q = demo_rvbpcamv
warning('This demo is deprecated. See folder fa instead')

n = 50;
m = 2;
d = 1;

randn('state', 6);
rand('state', 6);

mu = [0;0];%20 * randn(m,1);
W = [1; 0.5];%orth(randn(m,d)) * diag(10:-1:(10-d+1));
X = orth(randn(d,n)')' * sqrt(n);
s2 = 0.1;

% Generate normal data
Y = W*X + repmat(mu,1,n);

% Add noise
Yn = Y + sqrt(s2)*randn(m,n);

% Generate some outliers
Yno = Yn;
p = (rand(m,n) < 0.05);
%p = zeros(m,n) == 1;
%p(1,1) = true;
%p(2,2) = true;
%p(1,3) = true;
%p(:,4) = true;
Yno(p) = 3 + 2 * rand(sum(p(:)),1);

% Generate some missing values
Ynom = Yno;
pmv = (rand(m,n) < 0.0);
Ynom(pmv) = NaN;

% Use STANDARD PCA for non-outlier data
Q = vbpcamv(Yn, d);
W_pca = Q.W
X_pca = Q.X
mu_pca = Q.mu
s2_pca = 1/Q.tau
%[W_pca, X_pca, mu_pca, s2_pca] = pca_full(Yn, d);
Y_pca = W_pca*X_pca + repmat(mu_pca, 1, n);

% Run different algorithms
Q = vbpcamv(Ynom, d);
W_p = Q.W
X_p = Q.X
mu_p = Q.mu
s2_p = 1/Q.tau
%[W_p, X_p, mu_p, s2_p] = pca_full(Ynom, d);
[W_rp,X_rp,Sv_rp,mu_rp,nu_rp,s2_rp,U_rp] = rppcamv(Ynom, d, 'maxiters', 100);
init.tau = 1e2;
init.nu = 1;
init.mu = 0;%mean(Yno,2)
prior = [];
% $$$ prior.atau = 1e3;
% $$$ prior.atau = 1e0
prior.aw = 1e-10;
prior.bw = 1e-3;
results_rvb = rvbpcamv(Ynom, d, 'maxiters', 200, 'init', init, 'prior', ...
                       prior, 'startupdatehyper', 10, 'startrotate', 1, ...
                       'rotate', true);
W_rvb = results_rvb.W;
X_rvb = results_rvb.X;
S_rvb = results_rvb.Sv;
mu_rvb = results_rvb.mu;
nu_rvb = results_rvb.nu;

% Reconstruct
Y_p = W_p*X_p + repmat(mu_p,1,n);
Y_rp = W_rp*X_rp + repmat(mu_rp,1,n);
Y_rvb = W_rvb*X_rvb + repmat(mu_rvb,1,n);

% Error of the principal subspace
err_W_p = 180 * subspace(W_p, W_pca) / pi
err_W_rp = 180 * subspace(W_rp, W_pca) / pi
err_W_rvb = 180 * subspace(W_rvb, W_pca) / pi

plot_results(Ynom, []);
plot_results(Ynom, Y_p)
plot_results(Ynom, Y_rp)
plot_results(Ynom, Y_rvb)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
function plot_results(Y, Yh)
%subspace2d(Yh, [], Y)
figure
plot(Y(1,:), Y(2,:), 'x', 'markersize', 5, 'color', [.2 .2 .2]);
if ~isempty(Yh)
  hold on
  plot(Yh(1,:), Yh(2,:), 'ko', 'markersize', 3);
  nans = nan*ones(1,size(Y,2));
  plot([Y(1,:);Yh(1,:);nans], [Y(2,:);Yh(2,:);nans], ':', 'color', [.0 .0 .0]);
end
set(gca, 'DataAspectRatioMode', 'manual');
set(gca, 'DataAspectRatio', [1 1 1]);
set(gca, 'XTick', [], 'YTick', []);
set(gcf, 'units', 'centimeters');
pos = get(gcf, 'position');
set(gcf, 'position', [pos(1), pos(2), 5 5])
set(gcf, 'PaperPositionMode', 'auto', 'paperunits', 'centimeters', 'PaperSize', [5 5]);

xmg = 0.05 * (max(Y(1,:)) - min(Y(1,:)));
ymg = 0.05 * (max(Y(2,:)) - min(Y(2,:)));
xl = [min(Y(1,:))-xmg, max(Y(1,:))+xmg];
yl = [min(Y(2,:))-ymg, max(Y(2,:))+ymg];
set(gca, 'xlim', xl, 'ylim', yl);

% $$$ % Mark outliers
% $$$ hold on
% $$$ indeces = sum(p,1)>0;
% $$$ plot(VYno(1,indeces), VYno(2, indeces), 'go', 'MarkerSize', 8);

% $$$ % Mark observations with missing values
% $$$ indeces = sum(pmv,1)>0;
% $$$ scatter(VYno(1,indeces), VYno(2, indeces), 'yo');

return

