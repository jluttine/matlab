function [mapfigh, tsfigh] = mohsst5_plotexperiment(data, Q, isGP, comps)
% mohsst5_plotexperiment(data, Q, isGP, comps)

if nargin < 3
  isGP = true;
end

[LON,LAT] = meshgrid(data.longitude, data.latitude);
% $$$ [LAT,LON] = meshgrid(data.latitude, data.longitude);

variance_maps = false;

if isGP
  D = cols(Q.W);
else
  D = cols(Q.A);
end

if nargin < 4
  comps = 1:D;
end

% Rotate to PCA
weights = sqrt(cosd(data.coordinates(2,:)));
if isGP
  lands = colsum(~isnan(data.observations)) == 0;
  if rows(Q.W) == 1727 % t
    disp('The results are not ready? Fill the lands and use inverse weights..')
    Wall = zeros([rows(data.observations),D]);
    varWall = Wall;
    Wall(~lands,:) = Q.W;
    varWall(~lands,:) = Q.varW;
    Wall = bsxfun(@rdivide, Wall, weights(:));
    varWall = bsxfun(@rdivide, varWall, weights(:).^2);
  else
    Wall = Q.W;
    varWall = Q.varW;
  end
  % Remove land
  Wall(lands,:) = 0;
  varWall(lands,:) = 0;
  % Rotate using weights
  [Wall,varWall,Xall,varXall] = gprotate2pca(Wall,varWall,Q.X,Q.varX, ...
                                             weights(:).^2);
end


% Plot spatial patterns
mapproj('global-ellipse');
mapfigh = [];
for d=1:length(comps)
  if mod(d-1,8) == 0
    figure
    mapfigh = [mapfigh, gcf];
  end
  if variance_maps
    fig = 2*mod(d-1,8)+1;
    subplot(min(length(comps),8),2,fig);
  else
    fig = mod(d-1,8) + 1;
    subplot(1,min(length(comps),8),fig);
  end
  if isGP
    if numel(Wall(:,d)) ~= numel(LAT(:))
      error('Inputs do not correspond to interpolation points!')
      %[W,varW] = gppred(Q.pseudoW{d}, Q.Wp{d}, Q.CovWp{d}, [LON(:)'; ...
      %                    LAT(:)'], Q.logthetaW{d}, Q.covfuncW{d});
    else
      disp('Assuming that the inputs correspond to the interpolation points');
      W = Wall(:,comps(d));%Q.W(:,d);
      varW = varWall(:,comps(d));%Q.varW(:,d);
    end
  else
    W = Q.A(:,comps(d));
    varW = [];
  end
  set(gca, 'xtick', [], 'ytick', []);
  mappcolor(LON, LAT, reshape(W, size(LON)));
  h_cb = colorbar('SouthOutside');
  mapcoast;
  if ~isempty(varW) && variance_maps
    subplot(min(length(comps),8),2,fig+1);
    mappcolor(LON, LAT, reshape(2*sqrt(varW), size(LON)));
    colorbar;
    mapcoast;
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
  set( gca, 'Position', pos_ax )
  
  cb_pos = pos_ax;
  cb_pos(2) = 0.13;
  cb_pos(4) = 0.05;
  set(h_cb, 'Position', cb_pos);
  
  
end

set(gcf, 'Units', 'centimeters');
pos = get(gcf, 'Position');
figw = 35;
figh = 6;
set(gcf, 'Position', [pos(1) pos(2) figw figh]);
set(gcf, 'PaperPositionMode', 'auto', 'PaperSize', [figw figh]);


% Plot temporal features
input = data.time;
if isGP
  X = Xall;%Q.X;
  eX = 2*sqrt(varXall);2*sqrt(Q.varX);
else
  X = Q.S;
  eX = zeros(size(X));
  for j=1:cols(eX)
    eX(:,j) = 2*sqrt(diag(Q.Sv{j}));
  end
end
%if D > 8
years = 1800:2100; dates = [years(:), ones(length(years), 2)];
tsfigh = [];
for i=1:ceil((length(comps)-1)/8)
  i0 = (i-1)*8 + 1;
  i1 = i0 + min(7, length(comps)-(i-1)*8-1);
  hax = tsgpplot(input, X(comps(i0:i1),:)', eX(comps(i0:i1),:)');
  for j = 1:length(hax)
    xlim(hax(j), [min(input) max(input)]);
    set( hax(j), 'YTick', [] )
    set( hax(j), 'YTickLabel', [] )
    datetick(hax(j), 'x', 10, 'keeplimits');
    if j ~= length(hax)
      set(hax(j), 'xticklabel', []);
    end
  end
  tsfigh = [tsfigh, gcf];
end

set(gcf, 'Units', 'centimeters');
pos = get(gcf, 'Position');
figw = 35;
figh = 8;
set(gcf, 'Position', [pos(1) pos(2) figw figh]);
set(gcf, 'PaperPositionMode', 'auto', 'PaperSize', [figw figh]);










%else
%  tsgpplot(input, X', eX')
%end

% $$$ figure
% $$$ for d=1:D
% $$$   subplot(D,1,d);
% $$$   [W,varW] = gppred(Q_gp.pseudoW{d}, Q_gp.Wp{d}, Q_gp.CovWp{d}, [LON(:)'; ...
% $$$                       LAT(:)'], Q_gp.logthetaW{d}, Q_gp.covfuncW{d});
% $$$   m_pcolor(LON, LAT, reshape(W, size(LON)));
% $$$   shading('flat');
% $$$   mapcoast;
% $$$   cl = get(gca, 'clim');
% $$$   cl = max(abs(cl));
% $$$   set(gca, 'clim', [-cl cl]);
% $$$   mapcolormap;
% $$$   colorbar;
% $$$ end
