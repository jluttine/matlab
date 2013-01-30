function X = mohsst5_remove_land(X)

data = mohsst5_loaddata();

land = sum(~isnan(data.observations),2) == 0;

X(land,:) = [];