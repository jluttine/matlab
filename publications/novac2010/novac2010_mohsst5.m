function [data,Qgp,Qvb,ax] = novac2010_mohsst5(data, Qgp, Qvb)

if nargin < 1 || isempty(data)
  data = mohsst5_loaddata;
end

if nargin < 2 || isempty(Qgp)
  % GPPCA dates, e.g., 20091001 (crashed), 20090824, 20090820 (crashed),
  % 20090819 (crashed)
  Qgp = load(['mohsst5_weighted_testsize=20_D=80_Wcov=se_Wden=30_Xcov=se-' ...
              'per-cs_Xden=5_iter=200_date=20090824'], 'W', 'varW', 'X', ...
             'varX', 'tau', 'pseudoX', 'Xp', 'CovXp', 'logthetaX', ...
             'covfuncX');
% $$$   Qgp = load(['mohsst5_weighted_testsize=20_D=80_Wcov=se_Wden=30_Xcov=se-' ...
% $$$               'per-cs_Xden=5_iter=50_date=20090604'], 'W', 'varW', 'X', ...
% $$$              'varX', 'tau', 'loglikelihood');
end

if nargin < 3
  % VBPCA
  Qvb = load(['mohsst5_weighted_VBPCA_testsize=20_D=80_iter=50_date=' ...
              '20090601']);
  % Do rotation for some components
  Qvb.A(:,4) = -Qvb.A(:,4);
  Qvb.S(4,:) = -Qvb.S(4,:);

end


weights = sqrt(cosd(data.coordinates(2,:)));
lands = colsum(~isnan(data.observations)) == 0;
if rows(Qgp.W) == 1727 % t
  D = cols(Qgp.W);
  disp('The results are not ready? Fill the lands and use inverse weights..')
  W = zeros([rows(data.observations),D]);
  varW = W;
  W(~lands,:) = Qgp.W;
  varW(~lands,:) = Qgp.varW;
  Qgp.W = bsxfun(@rdivide, W, weights(:));
  Qgp.varW = bsxfun(@rdivide, varW, weights(:).^2);
end

% $$$ totally_missing_columns = sum( sum(~isnan(data.observations),1) == 0)
% $$$ return

[M,N] = size(data.observations);

% Choose what to show
plot_comps = false;
mixing_weights = true;

% $$$ 
% $$$ % Plot spatial and temporal components
% $$$ if plot_comps
% $$$   
% $$$   % GPFA
% $$$   [mapfigh, tsfigh] = mohsst5_plotexperiment(data, Qgp, true, 1:4);
% $$$ 
% $$$   % Saturate the second component
% $$$   warning('Saturating the second spatial component of GPFA!')
% $$$   ax = findobj(mapfigh, 'type', 'axes')
% $$$   set(ax(6), 'clim', [-0.7 0.7]);
% $$$   set(ax(4), 'clim', [-0.9 0.9]);
% $$$ %  return
% $$$   
% $$$ % $$$   print(mapfigh, '-depsc2', ['/home/jluttine/papers/2009NIPS/' ...
% $$$ % $$$                     'fig_experiment_loadings_gppca.eps']);
% $$$ % $$$   print(tsfigh, '-depsc2', ['/home/jluttine/papers/2009NIPS/' ...
% $$$ % $$$                     'fig_experiment_timeseries_gppca.eps']);
% $$$ 
% $$$   % VBPCA
% $$$   [mapfigh, tsfigh] = mohsst5_plotexperiment(data, Qvb, false, 1:4);
% $$$ 
% $$$   % Saturate the VBPCA components
% $$$   warning('Saturating VBPCA components');
% $$$   ax = findobj(mapfigh, 'type', 'axes')
% $$$   set(ax(8), 'clim', [-1 1]);
% $$$   set(ax(6), 'clim', [-0.8 0.8]);
% $$$   set(ax(4), 'clim', [-0.7 0.7]);
% $$$   set(ax(2), 'clim', [-0.7 0.7]);
% $$$ 
% $$$ % $$$   print(mapfigh, '-depsc2', ['/home/jluttine/papers/2009NIPS/' ...
% $$$ % $$$                     'fig_experiment_loadings_vbpca.eps']);
% $$$ % $$$   print(tsfigh, '-depsc2', ['/home/jluttine/papers/2009NIPS/' ...
% $$$ % $$$                     'fig_experiment_timeseries_vbpca.eps']);
% $$$ end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Mixing weights from the GP components to PCA components
if mixing_weights
  
  indsX = 1:N;

  %%%% FOR VBPCA %%%%
  
  if true
    
    % DON'T ROTATE VBPCA RESULTS BECAUSE THE COVARIANCE MATRIX Av IS
    % CORRUPTED!!!!
    % 
    % An adhoc solution:
    Av = covcell_to_covarray(Qvb.Av);
    Av(:) = 0;
    
    plot_pca_components(Qvb.A, Av, Qvb.S, covcell_to_covarray(Qvb.Sv), ...
                        data,'VBPCA', indsX, false);
    
  end

  %%%% FOR GPPCA %%%%

  %% Preprocess
  
  spatial_pca = false;
  
  if ~spatial_pca
    % First, remove bias and scale the GP components to unit variance
    X0 = bsxfun(@minus, Qgp.X, mean(Qgp.X,2));
    % Scale to unit variance
    xx = var(X0,1,2) + mean(Qgp.varX,2); % variance whitening
    sqrtxx = sqrt(xx);
    varX = bsxfun(@rdivide, Qgp.varX, xx);
    X = bsxfun(@rdivide, X0, sqrtxx);
    varW = bsxfun(@times, Qgp.varW, xx');
    W = bsxfun(@times, Qgp.W, sqrtxx');
  else
    disp('Plotting spatial PCA, is this correct??')
    % First, remove bias and scale the GP components to unit variance
    W0 = bsxfun(@minus, Qgp.W, mean(Qgp.W,1));
    % Scale to unit variance
    xx = mean(Qgp.X.^2,2) + mean(Qgp.varX,2); % 2nd moment whitening
    % xx = var(Qgp.X,1,2) + mean(Qgp.varX,2); % variance whitening
    sqrtxx = sqrt(xx);
    varX = bsxfun(@rdivide, Qgp.varX, xx);
    X = bsxfun(@rdivide, Qgp.X, sqrtxx);
    varW = bsxfun(@times, Qgp.varW, xx');
    W = bsxfun(@times, W0, sqrtxx');
  end
  
  % Variances to covariance matrices
  [M,D] = size(W);
  [D,N] = size(X);
  CovX = zeros([D D N]);
  for n=1:N
    CovX(:,:,n) = diag(varX(:,n));
  end
  CovW = zeros([D D M]);
  for m=1:M
    CovW(:,:,m) = diag(varW(m,:));
  end

  plot_pca_components(W,CovW,X, CovX,data,'GPPCA', indsX, true);


end



function [W,CovW,X,CovX,R,scaledR] = plot_pca_components(W,CovW,X,CovX, ...
                                                  data,method, indsX, rotate)

[M,D] = size(W);
[D,N] = size(X);

if nargin < 7 || isempty(indsX)
  indsX = 1:N;
end
if nargin < 8 || isempty(rotate)
  rotate = true;
end

% Ignore the not requested indices
X(:,setdiff(1:N,indsX)) = 0;
CovX(:,:,setdiff(1:N,indsX)) = 0;

% Ignore land areas
lands = colsum(~isnan(data.observations)) == 0;
W(lands,:) = 0;
CovW(:,:,lands) = 0;

weights2 = cosd(data.coordinates(2,:));
if rotate
  % Rotate to pca
  [W,CovW,X,CovX,R] = rotate_to_pca(W,CovW,X,CovX,weights2);
else
  R = eye(D);
end

% Evaluate explained variance for each component
WW = W'*diag(weights2)*W;
for m=1:M
  WW = WW + weights2(m) * CovW(:,:,m);
end

% Debug stuff:
%XX = X*X' + sum(CovX,3)
%error = mean(mean( abs(W*X - W_old*X_old) ))

%WW10 = WW(1:20,1:20)
S = sqrt(diag(diag(WW/size(W,1))));
scaledR = S*R;

if strcmpi(method, 'gppca')
  % Change the sign of some components
  R = eye(D);
  R(1,1) = -1;
  R(2,2) = -1;
  W = W * R;
  X = R * X;
end
if strcmpi(method, 'vbpca')
  % Change the sign of some components
  R = eye(D);
  R(3,3) = -1;
  W = W * R;
  X = R * X;
end

if true
  % Plot temporal PCA components
  Xpca = X;
  varXpca = zeros(size(X));
  for n=1:N
    varXpca(:,n) = diag(CovX(:,:,n));
  end
  eX = 2*sqrt(varXpca);
  hax = tsgpplot(data.time, Xpca(1:4,:)', eX(1:4,:)');

  % Publication style
  for j = 1:length(hax)
    xlim(hax(j), [min(data.time) max(data.time)]);
    set( hax(j), 'YTick', [] )
    set( hax(j), 'YTickLabel', [] )
    datetick(hax(j), 'x', 10, 'keeplimits');
    if j ~= length(hax)
      set(hax(j), 'xticklabel', []);
    end
  end
% $$$   set(gcf, 'Units', 'centimeters');
% $$$   pos = get(gcf, 'Position');
% $$$   figw = 35;
% $$$   figh = 8;
% $$$   set(gcf, 'Position', [pos(1) pos(2) figw figh]);
% $$$   set(gcf, 'PaperPositionMode', 'auto', 'PaperSize', [figw figh]);
  set_subplot_positions(hax, 4,1, [0.01 0.01 0.01 0.3], [0.01 0.01]);
  set_ticks_fontsize(hax, 5);
  set_figure_size(gcf, 13, 7);
  filename = sprintf(['/home/jluttine/thesis/figures/' ...
                      'novac2010_experiment_timeseries_%s.eps'], lower(method));
  print(gcf, '-depsc2', filename);
end

if true
  % Plot spatial PCA components
  figure
  mapproj('global-ellipse');
  [LON,LAT] = meshgrid(data.longitude, data.latitude);
  comps = 1:4;
  hax = [];
  h_cb = [];
  for d=comps
    subplot(1,4,d)
    set(gca, 'xtick', [], 'ytick', []);
    mappcolor(LON, LAT, reshape(W(:,d), size(LON)));
    h_cb(d) = colorbar('SouthOutside');
    hax(d) = gca;
    mapcoast;

    if strcmpi(method, 'vbpca')
      switch d
       case 1
        disp('Saturating vbpca 1-component')
        set(gca, 'clim', [-0.8 0.8]);
       case 2
        disp('Saturating vbpca 2-component')
        set(gca, 'clim', [-0.8 0.8]);
       case 4
        disp('Saturating vbpca 4-component')
        set(gca, 'clim', [-0.8 0.8]);
      end
    end
  end

  % Publication style:
  set_subplot_positions(hax, 1,4, [0.01 0.01 0.01 0.3], [0.01 0.01]);
  set_colorbar_position(h_cb, hax, 0.04, 0.04);
  set(h_cb, 'FontSize', 5);
  set_figure_size(gcf, 13, 2.3);
  filename = sprintf(['/home/jluttine/thesis/figures/' ...
                      'novac2010_experiment_loadings_%s.eps'], lower(method));
  print(gcf, '-depsc2', filename);

end