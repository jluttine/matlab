function X = mohsst5_remove_land(X)

data = mohsst5_loaddata();

land = sum(~isnan(data.observations),2) == 0;

Y = nan(length(land),size(X,2));
Y(~land,:) = X;
