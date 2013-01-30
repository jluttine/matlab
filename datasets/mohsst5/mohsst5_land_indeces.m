function ind = mohsst5_land_indeces()

data = mohsst5_loaddata();

land = sum(~isnan(data.observations),2) == 0;

vec = 1:length(land);
ind = vec(land);