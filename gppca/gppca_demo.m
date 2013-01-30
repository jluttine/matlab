function gppca_demo

% Good seeds: 13
randn('state', 16)
rand('state', 16)

inX = 100*(1:200);
inW = 10*rand(2,30);

D = 0;

% Trend
D = D + 1;
covfuncX{D} = @gpcov;
logthetaX{D} = log(365*(10+10*rand)); % 10-20 years
initthetaX{D} = logthetaX{D} - 0.1;%log(365*(10+10*rand));
densityX(D) = 1;

% Periodical
D = D + 1;
covfuncX{D} = {@gpcovProduct, @gpcov, @gpcovPeriodic};
%logthetaX{D} = log([5*365, 2, 4*365]); % 1 year
logthetaX{D} = log([100*365, 1.6, 3*365]); % decay 100 years, period 3 year
initthetaX{D} = logthetaX{D};%log([300*365, 2, 365*(3.5+1*rand)]);
densityX(D) = 1;

% Short-scale
D = D + 1;
covfuncX{D} = @gpcovPP;
logthetaX{D} = log(30*(24)); % 2 year
initthetaX{D} = logthetaX{D} - 0.1;%log(30*(20+20*rand));
densityX(D) = 1;

% Short-scale, almosta noise
D = D + 1;
covfuncX{D} = @gpcovPP;
logthetaX{D} = log(30*4); % 4 months
initthetaX{D} = logthetaX{D} - 0.1;%log(30*(3+6*rand));
densityX(D) = 1;

% Generate latent signals
for d=1:D
  X(d,:) = gprnd(inX, logthetaX{d}, covfuncX{d});
end

tsplot(X)
% $$$ return

scales = [5 8 10 8];
lengthscales = [5 2 2 0.3];
bias =   [3 4 0 0];
for d=1:D
% $$$   covfuncW{d} = {@(logtheta,x1,x2) gpcov(logtheta,x1,x2,@sqdistEuclidean)};
% $$$   logthetaW{d} = [log(2)];
% $$$   initthetaW{d} = [log(6)];
  covfuncW{d} = {@gpcovScale, @(logtheta,x1,x2) gpcov(logtheta,x1,x2,@ ...
                                                    sqdistEuclidean)};
  logthetaW{d} = [log(scales(d)); log(lengthscales(d))];
  initthetaW{d} = logthetaW{d} + 0.1;%[log(5); log(6)];
  densityW(d) = 1;
  W(:,d) = bias(d) + gprnd(inW, logthetaW{d}, covfuncW{d});
end

% $$$ size_inW = size(inW)
% $$$ size_W = size(W)

% Data
Y = W*X;

% Noise
Yn = Y + 1*randn(size(Y));

% Missing values
Ynm = Yn;
Itrain = rand(size(Yn))<0.7;
Ynm(Itrain) = nan;
Ytest = Yn;
Ytest(~Itrain) = nan;
Ytest(rand(size(Ytest))<0.87) = nan;

% Gap
Ygap = nan*Ynm;
Ynm(:,100:150) = nan;

% $$$ size_Y = size(Ynm)
% $$$ num_of_observations = sum(~isnan(Ynm(:)))
% $$$ 
% $$$ first_Y = Y(:,1:5)

Qgp = vbgppcamv_full(Ynm, D,...
                     inW, inX,...
                     covfuncW, initthetaW,...
                     covfuncX,initthetaX, ...
                     'maxiter',10, ...
                     'pseudodensityx', densityX, ...
                     'pseudodensityw', densityW, ...
                     'loglikelihood', true, ...
                     'updatehyper', [1 6 21 51 101 201 501], ...
                     'updatepseudox', false, ...
                     'updatepseudow', false, ...
                     'maxsearchx', 3, ...
                     'maxsearchw', 3, ...
                     'checkgradw', false, ...
                     'checkgradx', false);



Ygp = Qgp.W * Qgp.X;
varYgp = Qgp.W.^2*Qgp.varX + Qgp.varW*Qgp.X.^2 + Qgp.varW*Qgp.varX + 1/Qgp.tau; % CORRECT???
%trainrmse_gppca = rmse(Yobsw, Yh_gppca)
%testrmse_gppca = rmse(Ytestw, Yh_gppca)

% Plot latent signals
size_inX = size(inX)
size_X = size(Qgp.X)
size_varX = size(Qgp.varX)
hax = tsgpplot(inX', Qgp.X', 2*sqrt(Qgp.varX)');
for i=1:length(hax)
  %set(hax(i), 'xtick', [], 'ytick', []);
  axes(hax(i));
  line(100*[100, 100], [-1000 1000], 'Color', [0 0 0])
  line(100*[150, 150], [-1000 1000], 'Color', [0 0 0])
end


%print('-depsc2', '/home/jluttine/papers/2009NIPS/fig_artificial_latent.eps');

% Compare latent signals
% $$$ Qgp.X = bsxfun(@minus, Qgp.X, mean(Qgp.X,2));
% $$$ Qgp.varX = bsxfun(@rdivide, Qgp.varX, std(Qgp.X,1,2).^2);
% $$$ Qgp.X = bsxfun(@rdivide, Qgp.X, std(Qgp.X,1,2));
% $$$ X = bsxfun(@minus, X, mean(X,2));
% $$$ X = bsxfun(@rdivide, X, std(X,1,2));
% $$$ tsgpplot(inX', Qgp.X', 2*sqrt(Qgp.varX)');
% $$$ addtsplot(inX, X, 'r')

% Show predictive distribution
% $$$ hax = tsgpplot(inX', Ygp([1],:)', 2*sqrt(varYgp(1,:))');
% $$$ addtsplot(inX', Ynm([1],:), 'k+', 'MarkerSize', 10);
% $$$ addtsplot(inX', Ytest([1],:), 'ko', 'MarkerSize', 10);
M = size(Qgp.W,1);
random = randperm(M);
examples = random(1:4);
hax = tsgpplot(inX', Ygp(examples,:)', 2*sqrt(varYgp(examples,:))');
addtsplot(inX', Ynm(examples,:), 'k+', 'MarkerSize', 10);
addtsplot(inX', Ytest(examples,:), 'ko', 'MarkerSize', 10);
for i=1:length(hax)
  %set(hax(i), 'xtick', [], 'ytick', []);
  axes(hax(i));
  line(100*[100, 100], [-1000 1000], 'Color', [0 0 0])
  line(100*[150, 150], [-1000 1000], 'Color', [0 0 0])
end
%print('-depsc2', '/home/jluttine/papers/2009NIPS/fig_artificial_predictive.eps');

% Plot spatial components
[coord1,coord2] = meshgrid(linspace(0,10,50), linspace(0,10,50));
inWh = [coord1(:)'; coord2(:)'];
Mh = size(inWh,2);
Wh = zeros(Mh,D);
varWh = zeros(Mh,D);
figure
for d=1:D
  first = (d-1)*M + 1;
  last = first + M - 1;
  inds = first:last;
  [Wh(:,d), varWh(:,d)] = gppred(Qgp.inW, Qgp.W(:,d), Qgp.CovW(inds,inds), ...
                                 inWh, Qgp.logthetaW{d}, Qgp.covfuncW{d});
  
  % Plot mean map
  subplot(D,2,2*d-1)
  contourf(coord1,coord2,reshape(Wh(:,d), size(coord1)),20);
  hold on
  plot(inW(1,examples),inW(2,examples), 'kx', 'MarkerSize', 10);
  shading('flat');
  cl = get(gca, 'clim');
  cl = max(abs(cl));
  set(gca, 'clim', [-cl cl]);
  mapcolormap
  % Plot uncertainty map
  subplot(D,2,2*d);
  contourf(coord1,coord2,reshape(sqrt(varWh(:,d)), size(coord1)),20);
  hold on
  plot(inW(1,examples),inW(2,examples), 'kx', 'MarkerSize', 10);
  shading('flat');
  cl = get(gca, 'clim');
  cl = max(abs(cl));
  set(gca, 'clim', [-cl cl]);
  mapcolormap
end


% $$$ exp(Qgp.logthetaX{1})
% $$$ exp(Qgp.logthetaX{2})
% $$$ exp(Qgp.logthetaX{3})
