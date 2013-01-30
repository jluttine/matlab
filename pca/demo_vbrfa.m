% DEMO_VBRFA - A simple comparison of VB PCA and two robust extensions.
%
% The algorithms do not use ARD prior for the loadings W in order to keep
% the visualizations of the subspaces more clear. In addition, it helps
% avoiding pruning out relevant components too easily.
%
% NOTE: The robust algorithms might be less sensitive to the initialization
% if they estimate many components. With only one estimated component, they
% might kill all the components.
%
% This demo doesn't always show good results for robust algorithms. You may
% want to run a few times. In addition, sometimes there are several good
% interpretations of the data and the robust algorithms pick one.

% Last modified 2010-06-10
% Copyright (c) Jaakko Luttinen (jaakko.luttinen@tkk.fi)

function Q = demo_vbrfa(seed)

n = 50;
m = 2;
d = 1;
dh = d + 0; % number of components to estimate

if nargin >= 1
  randn('state', seed);
  rand('state', seed);
end

mu = [10;20];
W = [2; 1];
X = orth(randn(d,n)')' * sqrt(n);
s2 = 0.1;

% Generate normal data
Y = W*X + repmat(mu,1,n);

% Add noise
Yn = Y + sqrt(s2)*randn(m,n);

% Generate some outliers
Yno = Yn;
p = (rand(m,n) < 0.05);
Yno(p) = Yno(p) + (5 + 8 * rand(sum(p(:)),1));

% Generate some missing values
Ynom = Yno;
pmv = (rand(m,n) < 0.0);
Ynom(pmv) = NaN;

% Run different algorithms
options.init.tau = 1e1*ones(m,1); % the initialization of this can be crucial..
options.init.nu = 1;
options.prior.a_alpha = 1e10; % "fix" w to non-informative for W by
options.prior.b_alpha = 1e15; % setting strong prior for w
options.update_nu = 1;
options.update_alpha = 1;
options.update_beta = 1;
options.rotate = false;
options.common_nu = false;
options.common_tau = true;
options.maxiter = 100;

% PCA
disp('Run PCA')
options.robustness = 'none';
results_p = vbrfa(Ynom, dh, options);
W_p = results_p.W;
X_p = results_p.X;
Mu_p = results_p.Mu;
nu_p = results_p.nu;

% Robust PCA with multivariate Student-t
disp('Run robust PCA with multivariate Student-t')
options.robustness = 'multivariate-t';
results_rp = vbrfa(Ynom, dh, options);
W_rp = results_rp.W;
X_rp = results_rp.X;
Mu_rp = results_rp.Mu;
nu_rp = results_rp.nu;

options.maxiter = 1000;

% Robust PCA with independent Student-t
disp('Run robust PCA with independent Student-t')
options.robustness = 'independent-t';
results_rvb = vbrfa(Ynom, dh, options);
W_rvb = results_rvb.W;
X_rvb = results_rvb.X;
Mu_rvb = results_rvb.Mu;
nu_rvb = results_rvb.nu;

% Reconstruct
Y_p = bsxfun(@plus, W_p*X_p, Mu_p);
Y_rp = bsxfun(@plus, W_rp*X_rp, Mu_rp);
Y_rvb = bsxfun(@plus, W_rvb*X_rvb, Mu_rvb);

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
% $$$ pos = get(gcf, 'position');
% $$$ set(gcf, 'position', [pos(1), pos(2), 5 5])
% $$$ set(gcf, 'PaperPositionMode', 'auto', 'paperunits', 'centimeters', 'PaperSize', [5 5]);

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

