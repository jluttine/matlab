function b = mohsst5_is_land_index(ind)

data = mohsst5_loaddata();

land = sum(~isnan(data.observations),2) == 0;

b = land(ind);