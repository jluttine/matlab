function [data,Qgp,Qvb,ax] = mohsst5_paperplots(data, Qgp, Qvb)

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
rmse = false;
explained_variance = false;
comps = false;
ranked = false;
mixing_weights = true;
pred = false;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RMSE
if rmse
% $$$   rmse_gp = mohsst5_rmse_sparsepart(Ygp, data)
% $$$   rmse_vb = mohsst5_rmse_sparsepart(Yvb, data)

  % Test set
  Itest = load('ind20test.mat');
  Itest = Itest.Itest;
  testset = data.observations;
  testset(~Itest) = nan;

  % Weighted squared errors
  Ygp_se = bsxfun(@times, weights(:).^2, (testset - Qgp.W * Qgp.X) .^ 2);
  Yvb_se = bsxfun(@times, weights(:).^2, (testset - bsxfun(@plus, Qvb.A*Qvb.S, ...
                                                    Qvb.Mu)) .^ 2);
  Obs = ~isnan(testset);
  % Weight of each individual value (by number of observations)
  Obs_weighted = bsxfun(@times, weights(:).^2, Obs);
  Nh = sum(Obs_weighted, 2);
  Mh = sum(Obs_weighted, 1);
  NMh = sum(sum( Obs_weighted ));
  num_of_testvalues = sum( Obs(:) )
  
  % Zeros for missing and non-test values
  Ygp_se(~Obs) = 0;
  Yvb_se(~Obs) = 0;

  % Total RMSEs
  rmse_gp = sqrt(sum(Ygp_se(:))/NMh)
  rmse_vb = sqrt(sum(Yvb_se(:))/NMh)
  
  % RMSE temporally
  figure
  subplot(3,1,1)
  plot(M-sum(Obs,1))
  % RMSE as a function of time
  subplot(3,1,2)
  gp_temporal_rmse = sqrt(sum(Ygp_se,1)./Mh);
  plot(gp_temporal_rmse, 'k')
  hold on
  vb_temporal_rmse = sqrt(sum(Yvb_se,1)./Mh);
  plot(vb_temporal_rmse, 'r')
  % Smoother temporal curves:
  w = 48; % number of months to average over
  gp_smooth_temporal_mse = moving_average(Mh.*gp_temporal_rmse.^2, w/2) ./ ...
      moving_average(Mh,w/2);
  subplot(3,1,3)
  plot(sqrt(gp_smooth_temporal_mse), 'k');
  hold on
  vb_smooth_temporal_mse = moving_average(Mh.*vb_temporal_rmse.^2, w/2) ./ ...
      moving_average(Mh,w/2);
  subplot(3,1,3)
  plot(sqrt(vb_smooth_temporal_mse), 'r');

  % RMSE spatially
  figure
  mapproj('global-ellipse');
  [LON,LAT] = meshgrid(data.longitude, data.latitude);
  % Number of missing values
  subplot(2,2,1)
  set(gca, 'xtick', [], 'ytick', []);
  Nh(Nh==0) = eps;
  values = N-sum(Obs,2);
  values(lands) = nan;
  mappcolor(LON, LAT, reshape(values, size(LON)));
  set(gca, 'clim', [min(values) max(values)]);
  h_cb = colorbar('SouthOutside');
  mapcoast;
  % GP
  subplot(2,2,3)
  set(gca, 'xtick', [], 'ytick', []);
  Nh(Nh==0) = eps;
  values = sqrt(sum(Ygp_se,2)./Nh);
  values(lands) = nan;
  mappcolor(LON, LAT, reshape(values, size(LON)));
  set(gca, 'clim', [min(values) 1]);%max(values)]);
  h_cb = colorbar('SouthOutside');
  mapcoast;
  % PCA
  subplot(2,2,4)
  set(gca, 'xtick', [], 'ytick', []);
  Nh(Nh==0) = eps;
  values = sqrt(sum(Yvb_se,2)./Nh);
  values(lands) = nan;
  mappcolor(LON, LAT, reshape(values, size(LON)));
  set(gca, 'clim', [min(values) 1]);%max(values)]);
  h_cb = colorbar('SouthOutside');
  mapcoast;

  
end

%%%%%%%%%%%%%%
if explained_variance
% $$$   % VBPCA
% $$$   % WTF?!?!?!!! Qvb.Av is 1000*eye matrices?!?!?!!!
% $$$   X = bsxfun(@minus, Qvb.S, mean(Qvb.S,2)); % zero mean
% $$$   XX = X*X' + sum(covcell_to_covarray(Qvb.Sv),3);
% $$$   X2 = X*X';
% $$$   W = Qvb.A;
% $$$   WW = W'*W + sum(covcell_to_covarray(Qvb.Av),3);
% $$$   W2 = W'*W;
% $$$   Yh = W*X;
% $$$   %  varY = diag(XX*WW - X2*W2) / (N*M);
% $$$   meanY = mean(Yh(:));
% $$$   varY = diag(XX*WW) / (N*M) - mean(Yh(:))^2;
% $$$   figure
% $$$   plot(varY)
% $$$   varY_total_vbpca = sum(varY)
  
  
  % Ignore land areas
  W = Qgp.W;
  varW = Qgp.varW;
  lands = colsum(~isnan(data.observations)) == 0;
  W(lands,:) = 0;
  varW(lands,:) = 0;
  
  % Use weighted variance..
  weights2 = weights(:).^2;

  % "THE QUESTION": What is the amount of temporal variance each
  % component explains?
  
  %% JAAKKO'S WAY

  % GPFA (by zero meaning X this shows explained TEMPORAL variance)
  X0 = Qgp.X; disp('Not zero meaning X!!'); % do not zero mean
  %X0 = bsxfun(@minus, Qgp.X, mean(Qgp.X,2)); % zero mean
  % Scale to temporal components to unit variance
  xx = var(X0,1,2) + mean(Qgp.varX,2); % variance whitening
  sqrtxx = sqrt(xx);
  varX = bsxfun(@rdivide, Qgp.varX, xx);
  X = bsxfun(@rdivide, X0, sqrtxx);
  varW = bsxfun(@times, varW, xx');
  W = bsxfun(@times, W, sqrtxx');
  % Evaluate variances
  XX = X*X' + diag(sum(varX,2));
  WW = W'*diag(weights2)*W + diag(weights2'*varW);
  %Yh = W*X;
  mX = mean(X,2);
  varY = diag(XX*WW) / (N*M) - diag((W'*diag(weights2)*W)*(mX*mX')) / M;
%  varY = diag(XX*WW) / (N*M) - mean(Yh(:))^2 / length(XX);
  % Plot variance curve
  figure
  bar(varY)
  title('Explained variance by Jaakko')
  % Plot components explaining most variance
  varY_total_gpfa = sum(varY)
  [sorted_variances, inds] = sort(varY, 'descend');
  variances_per_components = [sorted_variances(:), inds(:)]
  tsgpplot(data.time(:), X(inds(1:8),:)', 2*sqrt(varX(inds(1:8),:))');

  %% ALEXANDER'S WAY
  
  X0 = bsxfun(@minus, Qgp.X, mean(Qgp.X,2)); % zero mean
  XX = X0*X0' + diag(sum(Qgp.varX,2));
  [V,D] = svd(XX/N);
  RX = V*D^(-0.5)*V'; % rotation
  X0_rotated = RX * X0;
  varX_rotated = RX * Qgp.varX;
  XX_rotated = X0_rotated*X0_rotated' + diag(sum(varX_rotated,2));
  RW = V*D^(0.5)*V';
  W_rotated = W * RW;
  varW_rotated = varW * RW;
  WW_rotated = W_rotated'*diag(weights2)*W_rotated + diag(weights2'*varW_rotated);
  variances = diag(WW_rotated) / M;
  figure
  bar(variances)
  title('Explained variance by Alexander')
  varianceY_total = sum(variances)
  
end

% Plot spatial and temporal components
if comps
  
  % GPFA
  [mapfigh, tsfigh] = mohsst5_plotexperiment(data, Qgp, true, 1:4);

  % Saturate the second component
  warning('Saturating the second spatial component of GPFA!')
  ax = findobj(mapfigh, 'type', 'axes')
  set(ax(6), 'clim', [-0.7 0.7]);
  set(ax(4), 'clim', [-0.9 0.9]);
%  return
  
% $$$   print(mapfigh, '-depsc2', ['/home/jluttine/papers/2009NIPS/' ...
% $$$                     'fig_experiment_loadings_gppca.eps']);
% $$$   print(tsfigh, '-depsc2', ['/home/jluttine/papers/2009NIPS/' ...
% $$$                     'fig_experiment_timeseries_gppca.eps']);

  % VBPCA
  [mapfigh, tsfigh] = mohsst5_plotexperiment(data, Qvb, false, 1:4);

  % Saturate the VBPCA components
  warning('Saturating VBPCA components');
  ax = findobj(mapfigh, 'type', 'axes')
  set(ax(8), 'clim', [-1 1]);
  set(ax(6), 'clim', [-0.8 0.8]);
  set(ax(4), 'clim', [-0.7 0.7]);
  set(ax(2), 'clim', [-0.7 0.7]);

% $$$   print(mapfigh, '-depsc2', ['/home/jluttine/papers/2009NIPS/' ...
% $$$                     'fig_experiment_loadings_vbpca.eps']);
% $$$   print(tsfigh, '-depsc2', ['/home/jluttine/papers/2009NIPS/' ...
% $$$                     'fig_experiment_timeseries_vbpca.eps']);
end

% Plot ranked GP components
if ranked

  error(['This is incorrect (Does not take into account correlations ' ...
         'between the components! Check out the explained variance part!'])
         
  % First, remove bias and scale the GP components to unit variance
  bias = mean(Qgp.X,2);
  X0 = bsxfun(@minus, Qgp.X, bias);
  % Scale to unit variance
  xx = var(X0,1,2) + mean(Qgp.varX,2); % variance whitening
  sqrtxx = sqrt(xx);
  varX = bsxfun(@rdivide, Qgp.varX, xx);
  X = bsxfun(@rdivide, X0, sqrtxx);
  varW = bsxfun(@times, Qgp.varW, xx');
  W = bsxfun(@times, Qgp.W, sqrtxx');
  
  % Ignore land areas
  lands = colsum(~isnan(data.observations)) == 0;
  W(lands,:) = 0;
  varW(lands,:) = 0;
  
  weights2 = cosd(data.coordinates(2,:));
  WW = W' * diag(weights2) * W + diag(weights2(:)'*varW);
  [vars,ind] = sort(diag(WW), 'descend');
  variances = vars(:)
  
  % Show indeces
  indeces = ind(1:8)
  Dh = length(indeces);
  
  if false
    % Predict a little bit?
    predtimes = data.time(end) + 30:30:900;
    Nh = length(predtimes);
    Xh = zeros(Dh, Nh);
    Qgp.CovXp
    size(Qgp.CovXp{1})
    for d=indeces
      % Predict
      [Xh(d,:), varXh(d,:)] = gppred(Qgp.pseudoX{d}, Qgp.Xp{d}, Qgp.CovXp{d}, ...
                                     predtimes, Qgp.logthetaX{d}, ...
                                     Qgp.covfuncX{d});
      % Transform as in preprocessing step
      Xh(d,:) = Xh(d,:) - bias(d);
      Xh(d,:) = Xh(d,:) / sqrtxx(d);
      varXh(d,:) = varXh(d,:) / xx(d);
    end

    plottimes = [data.time(:); predtimes(:)];
    Xplot = [X(indeces,:), Xh];
    varXplot = [varX(indeces,:), varXh];
  else
    plottimes = data.time(:);
    Xplot = X(indeces,:);
    varXplot = varX(indeces,:);
  end  
  
  % Plot temporal
  tsgpplot(plottimes, Xplot', 2*sqrt(varXplot)');
  
  % Plot spatial
  figure
  mapproj('global-ellipse');
  [LON,LAT] = meshgrid(data.longitude, data.latitude);
  for ind=1:Dh
    subplot(ceil(Dh/4),4,ind)
    set(gca, 'xtick', [], 'ytick', []);
    mappcolor(LON, LAT, reshape(W(:,indeces(ind)), size(LON)));
    h_cb = colorbar('SouthOutside');
    mapcoast;
  end
  
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Mixing weights from the GP components to PCA components
if mixing_weights
  
  if false % use only end part of the data
    indsX = 1300:N;
  else
    indsX = 1:N;
  end
  
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
  
  % First, remove bias and scale the GP components to unit variance
  X0 = bsxfun(@minus, Qgp.X, mean(Qgp.X,2));
  % Scale to unit variance
  xx = var(X0,1,2) + mean(Qgp.varX,2); % variance whitening
  sqrtxx = sqrt(xx);
  varX = bsxfun(@rdivide, Qgp.varX, xx);
  X = bsxfun(@rdivide, X0, sqrtxx);
  varW = bsxfun(@times, Qgp.varW, xx');
  W = bsxfun(@times, Qgp.W, sqrtxx');
  
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


% $$$   [Wpca,CovWpca,Xpca,CovXpca,R,scaledR] = plot_pca_components(W,CovW,X, ...
% $$$                                                     CovX,data,'GPPCA');

% $$$   % Rank components
% $$$   %[explained_var, ind] = sort(mean(scaledR.^2,1), 'descend')
% $$$   % Use another ranking
% $$$   % Ignore land areas
% $$$   lands = colsum(~isnan(data.observations)) == 0;
% $$$   W(lands,:) = 0;
% $$$   CovW(:,:,lands) = 0;
% $$$   weights2 = cosd(data.coordinates(2,:));
% $$$   [explained_var, ind] = sort(sqrt(weights2(:)'*(W.^2+varW)/M), 'descend')
% $$$   
% $$$   %compare = [explained_var(:), explained_var2(:), ind(:), ind2(:)]
% $$$   
% $$$   % Plot ranked GP components
% $$$   eX = 2*sqrt(varX);
% $$$   hax = tsgpplot(data.time, X(ind(1:4),:)', eX(ind(1:4),:)');
  
end

% Show predictions
if pred
  % Load test set indeces
  Itest = load('ind20test');
  Itest = Itest.Itest;

  % Test and train sets
  Ytest = data.observations;
  Ytest(~Itest) = nan;
  Ytrain = data.observations;
  Ytrain(Itest) = nan;

  % Reconstructions
  Ygp = Qgp.W * Qgp.X;
  varYgp = Qgp.W.^2*Qgp.varX + Qgp.varW*Qgp.X.^2 + Qgp.varW*Qgp.varX + 1/Qgp.tau;

  % Select stations
  inds = [300, 450, 820, 1200];

  Yvb = bsxfun(@plus, Qvb.A*Qvb.S, Qvb.Mu);
  varYvb = zeros(size(Yvb));
  for n=1:cols(Yvb)
    for m=inds
      varYvb(m,n) = Qvb.A(m,:)*Qvb.Sv{n}*Qvb.A(m,:)' + ...
          Qvb.S(:,n)'*Qvb.Av{m}*Qvb.S(:,n) + ...
          traceprod(Qvb.Sv{n}, Qvb.Av{m}) + ...
          Qvb.V;
    end
  end
  
  % Show predictive distribution for GP
  tsgpplot(data.time', Ygp(inds,:)', 2*sqrt(varYgp(inds,:))');
  %addtsplot(data.time, Ytrain(inds,:), 'r.');
  addtsplot(data.time, Ytest(inds,:), 'bx');

  % Show predictive distribution for VB
  tsgpplot(data.time', Yvb(inds,:)', 2*sqrt(varYvb(inds,:))');
  %addtsplot(data.time, Ytrain(inds,:), 'r.');
  addtsplot(data.time, Ytest(inds,:), 'bx');

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

if false
  % Plot mixing weights
  figure
  imagesc(abs(scaledR))
  xlabel([method, ' components'])
  ylabel('PCA components')
  title(['Mixing weights for ', method, ' (absolute values)']);
end

if false
  % Plot explained variances
  %variances = diag(WW) / M;
  %variances = diag(R'*WW*R)) / M;
  variances = diag(inv(R)*WW*R) / M;
  figure
  bar(variances)
  title('Explained variances of GPFA components in plot-pca-components')
  total_variance = sum(variances)
end

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
  set(gcf, 'Units', 'centimeters');
  pos = get(gcf, 'Position');
  figw = 35;
  figh = 8;
  set(gcf, 'Position', [pos(1) pos(2) figw figh]);
  set(gcf, 'PaperPositionMode', 'auto', 'PaperSize', [figw figh]);
  filename = sprintf(['/home/jluttine/papers/2009NIPS/' ...
                      'fig_experiment_timeseries_%s.eps'], lower(method));
  print(gcf, '-depsc2', filename);
end

if true
  % Plot spatial PCA components
  figure
  mapproj('global-ellipse');
  [LON,LAT] = meshgrid(data.longitude, data.latitude);
  comps = 1:4;
  for d=comps
    subplot(1,4,d)
    set(gca, 'xtick', [], 'ytick', []);
    mappcolor(LON, LAT, reshape(W(:,d), size(LON)));
    h_cb = colorbar('SouthOutside');
    mapcoast;

    % Publication style:
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
    n = length(comps);
    pos_ax = get( gca, 'Position' );
    hgt_ax = (0.95 - 0.1) / ( n + (n-1)*0.1 );
    hgt_sp = hgt_ax * 0.1;
    pos_ax(2) = 0.2;
    pos_ax(1) = d*( hgt_ax + hgt_sp ) - hgt_ax;
    %pos_ax(1) = 0.95 - (d-1)*( hgt_ax + hgt_sp ) - hgt_ax;
    pos_ax(4) = 0.74;
    pos_ax(3) = hgt_ax;
    set(gca, 'Color', 'none');
    set( gca, 'Position', pos_ax )
    cb_pos = pos_ax;
    cb_pos(2) = 0.13;
    cb_pos(4) = 0.05;
    set(h_cb, 'Position', cb_pos);
  end
  
  % Publication style
  set(gcf, 'Color', 'none');
  set(gcf, 'Units', 'centimeters');
  pos = get(gcf, 'Position');
  figw = 35;
  figh = 6;
  set(gcf, 'Position', [pos(1) pos(2) figw figh]);
  set(gcf, 'PaperPositionMode', 'auto', 'PaperSize', [figw figh]);
  filename = sprintf(['/home/jluttine/papers/2009NIPS/' ...
                      'fig_experiment_loadings_%s.eps'], lower(method));
  print(gcf, '-depsc2', filename);

end