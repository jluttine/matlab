% e = mohsst5_performance_rmsew(Yh, I)

function e = mohsst5_performance_rmsew(Yh, I)

load('mohsst5_data.mat', 'observations');

Y = observations;

if all(size(Yh) == [1727 1632])
  % The seas have been removed
  disp('Adding the land areas to the matrix.')
  sea = sum(~isnan(Y),2) > 0;
  F = nan(size(Y));
  F(sea,:) = Yh;
  Yh = F;
  Itmp = false(size(Y));
  Itmp(sea,:) = I;
  I = Itmp;
end

Y(~I) = nan;

e = mohsst5_rmsew(Y,Yh);
