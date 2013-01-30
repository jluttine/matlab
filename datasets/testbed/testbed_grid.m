% [LON,LAT] = TESTBED_GRID(RES_LON, RES_LAT)

% Copyright (c) 2010 Jaakko Luttinen

function [lon,lat] = testbed_grid(res_lon, res_lat)

[minlon,maxlon,minlat,maxlat] = testbed_projection();
lon = linspace(minlon, maxlon, res_lon);
lat = linspace(minlat, maxlat, res_lat);
