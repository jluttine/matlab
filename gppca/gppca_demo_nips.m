%function gppca_demo_nips


% Good seeds: 16, 20
seed = 16
randn('state', seed)
rand('state', seed)

inX = 100*(1:200);
inW = 10*rand(2,30);

N = size(inX,2);
M = size(inW,2);

D = 0;

% Trend
D = D + 1;
covfuncX{D} = @gpcov;
logthetaX{D} = log(365*(10+10*rand)); % 10-20 years
initthetaX{D} = logthetaX{D} + 1;%log(365*(10+10*rand));
densityX(D) = 1;

% Periodical (it is difficult to learn from a badly initialized
% periodicity :( )
D = D + 1;
covfuncX{D} = {@gpcovProduct, @gpcov, @gpcovPeriodic};
%logthetaX{D} = log([5*365, 2, 4*365]); % 1 year
logthetaX{D} = log([250*365, 1.6, 3*365]); % decay 100 years, period 3 year
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
X = zeros(D,N);
for d=1:D
  X(d,:) = gprnd(inX, logthetaX{d}, covfuncX{d});
end

scales = [2^2 3^2 3^2 1.0^2];
lengthscales = [5 1.5 1.5 0.9];
bias =   [3 0 0 0];
% Working set
% $$$ scales = [2^2 3^2 5^2 1.0^2];
% $$$ lengthscales = [5 2.5 1.5 0.9];
% $$$ bias =   [3 0 0 0];
W = zeros(M,D);
for d=1:D
% $$$   covfuncW{d} = {@(logtheta,x1,x2) gpcov(logtheta,x1,x2,@sqdistEuclidean)};
% $$$   logthetaW{d} = [log(2)];
% $$$   initthetaW{d} = [log(6)];
  covfuncW{d} = {@gpcovScale, @(logtheta,x1,x2) gpcov(logtheta,x1,x2,@ ...
                                                    sqdistEuclidean)};
  logthetaW{d} = [log(scales(d)); log(lengthscales(d))];
  initthetaW{d} = logthetaW{d} - 0.1;%[log(5); log(6)];
  densityW(d) = 1;
  W(:,d) = bias(d) + gprnd(inW, logthetaW{d}, covfuncW{d});
end


% $$$ size_inW = size(inW)
% $$$ size_W = size(W)

% Data
Y = W*X;
[M,N] = size(Y);

% Noise
Yn = Y + 1*randn(size(Y));

M = size(W,1);
random = randperm(M);
examples = random([4 9 14 15]);

% $$$ tsplot(X)
% $$$ tsplot(W')
% $$$ tsplot(Y(examples,:));
% $$$ return

% Missing values
Ynm = Yn;
Itrain = rand(size(Yn))<0.9;
Ynm(Itrain) = nan;
Ytest = Yn;
Ytest(~Itrain) = nan;
Ytest(rand(size(Ytest))<0.87) = nan;

% Create a gap
Ygap = nan*Ynm;
Ynm(:,100:150) = nan;

% $$$ size_Y = size(Ynm)
% $$$ num_of_observations = sum(~isnan(Ynm(:)))
% $$$ 
% $$$ first_Y = Y(:,1:5)

num_of_observations = sum(sum( ~isnan(Ynm) ))
%return

Qgp = vbgppcamv_full(Ynm, D,...
                     inW, inX,...
                     covfuncW, initthetaW,...
                     covfuncX,initthetaX, ...
                     'maxiter',50, ...
                     'pseudodensityx', densityX, ...
                     'pseudodensityw', densityW, ...
                     'loglikelihood', true, ...
                     'updatehyper', [3 5 10 20 50 80 100], ...
                     'updatepseudox', false, ...
                     'updatepseudow', false, ...
                     'maxsearchx', 3, ...
                     'maxsearchw', 3, ...
                     'checkgradw', false, ...
                     'checkgradx', false);

Ygp = Qgp.W * Qgp.X;
if false
  disp('Using inaccurate prediction variance.');
  varYgp = Qgp.W.^2*Qgp.varX + Qgp.varW*Qgp.X.^2 + Qgp.varW*Qgp.varX + ...
           1/Qgp.tau;
else
  disp('Using accurate prediction variance.')
  varYgp = zeros(size(Ygp));
  for m=1:M
    CovWm = Qgp.CovW(Qgp.indsW(m,:),Qgp.indsW(m,:));
    for n=1:N
      CovXn = Qgp.CovX(Qgp.indsX(:,n),Qgp.indsX(:,n));
      varYgp(m,n) = Qgp.W(m,:)*CovXn*Qgp.W(m,:)' + Qgp.X(:,n)'*CovWm*Qgp.X(:,n) ...
          + traceprod(CovXn,CovWm); 
    end
  end
  varYgp = varYgp + 1/Qgp.tau;
end

% Plot latent signals
hax = tsgpplot(inX', Qgp.X', 2*sqrt(Qgp.varX)');
for i=1:length(hax)
  set(hax(i), 'xtick', [], 'ytick', []);
  axes(hax(i));
  line(100*[100, 100], [-1000 1000], 'Color', [0 0 0])
  line(100*[150, 150], [-1000 1000], 'Color', [0 0 0])
  lab = sprintf('x_{%d}(t)', i);
  ylabel(lab);
end
xlabel('time, t');
set_figure_size(7, 6);
set_label_fontsize(hax, 6);
%print('-depsc2', '/home/jluttine/papers/2009NIPS/fig_artificial_latent');

% Compare latent signals
% $$$ Qgp.X = bsxfun(@minus, Qgp.X, mean(Qgp.X,2));
% $$$ Qgp.varX = bsxfun(@rdivide, Qgp.varX, std(Qgp.X,1,2).^2);
% $$$ Qgp.X = bsxfun(@rdivide, Qgp.X, std(Qgp.X,1,2));
% $$$ X = bsxfun(@minus, X, mean(X,2));
% $$$ X = bsxfun(@rdivide, X, std(X,1,2));
% $$$ tsgpplot(inX', Qgp.X', 2*sqrt(Qgp.varX)');
% $$$ addtsplot(inX, X, 'r')

% Show predictive distribution
hax = tsgpplot(inX', Ygp(examples,:)', 2*sqrt(varYgp(examples,:))');
addtsplot(inX, Ynm(examples,:), '+', 'MarkerSize', 4, 'Color', [0.8 ...
                    0 0]);
addtsplot(inX, Ytest(examples,:), 'o', 'MarkerSize', 4, 'Color', [0 ...
                    0 0.8]);
for i=1:length(hax)
  set(hax(i), 'xtick', [], 'ytick', []);
  axes(hax(i));
  line(100*[100, 100], [-1000 1000], 'Color', [0 0 0])
  line(100*[150, 150], [-1000 1000], 'Color', [0 0 0])
  lab = sprintf('y_{%d}(t)', examples(i));
  ylabel(lab);
end
xlabel('time, t');
set_figure_size(7, 6);
set_label_fontsize(hax, 6);
%print('-depsc2', ['/home/jluttine/papers/2009NIPS/' ...
%                    'fig_artificial_predictive']);

exp(Qgp.logthetaX{1})
exp(Qgp.logthetaX{2})
exp(Qgp.logthetaX{3})






% Plot spatial components
[coord1,coord2] = meshgrid(linspace(0,10,50), linspace(0,10,50));
inWh = [coord1(:)'; coord2(:)'];
Mh = size(inWh,2);
Wh = zeros(Mh,D);
varWh = zeros(Mh,D);
plot_spotlight = true;
figure
hax = [];
for d=1:D
  %first = (d-1)*M + 1;
  %last = first + M - 1;
  %inds = first:last;
  inds = Qgp.indsW(:,d);
  [Wh(:,d), varWh(:,d)] = gppred(Qgp.inW, Qgp.W(:,d), Qgp.CovW(inds,inds), ...
                                 inWh, Qgp.logthetaW{d}, Qgp.covfuncW{d});
  
  % Plot mean map
  if plot_spotlight
    subplot(2,D/2,d)
  else
    subplot(2,D,d)
  end
  contourf(coord1,coord2,reshape(Wh(:,d), size(coord1)),20);
  hold on
  plot(inW(1,:),inW(2,:), 'kx', 'MarkerSize', 8);
  %  plot(inW(1,examples),inW(2,examples), 'kx', 'MarkerSize', 10);
  shading('flat');
  cl = get(gca, 'clim');
  cl = max(abs(cl));
  set(gca, 'clim', [-cl cl]);
  map_colormap;
  set(gca, 'xtick', [], 'ytick', []);
  %pbaspect([1 1 1]);
  h_cb = colorbar('SouthOutside');
  if plot_spotlight
    cb_height = 0.03;
    cb_textheight = 0.07;
    cb_sep = 0.02;
    hax(d) = gca;
    set_subplot_positions(hax, 2, 2, [0.01 0.01 0.01 (cb_height+ ...
                                                      cb_textheight+cb_sep)], ...
                               [0.02 (cb_sep+cb_height+cb_textheight)]);
    pos_ax = get( gca, 'Position');
    cb_pos = pos_ax;
    cb_pos(2) = pos_ax(2) - cb_height - cb_sep;
    cb_pos(4) = cb_height;
    set(h_cb, 'Position', cb_pos);
  end
  
  if ~plot_spotlight
    % Plot uncertainty map
    subplot(2,D,D+d);
    contourf(coord1,coord2,reshape(sqrt(varWh(:,d)), size(coord1)),20);
    hold on
    plot(inW(1,:),inW(2,:), 'kx', 'MarkerSize', 8);
    %  plot(inW(1,examples),inW(2,examples), 'kx', 'MarkerSize', 10);
    shading('flat');
    cl = get(gca, 'clim');
    cl = max(abs(cl));
    set(gca, 'clim', [-cl cl]);
    map_colormap;
    set(gca, 'xtick', [], 'ytick', []);
    %pbaspect([1 1 1]);
    h_cb = colorbar('SouthOutside');
  end
  
end
if plot_spotlight
  set_figure_size(10, 11.5);
%  print(gcf, '-depsc2', ['/home/jluttine/papers/2009NIPS/' ...
%                      'fig_artificial_loadings']);
end