%
% [lon, lat] = lonlat2standard(lon, lat)
% [coord] = lonlat2standard(coord)
function [lon, lat] = lonlat2standard(lon, lat)

if nargin == 1
  lat = lon(2,:);
  lon = lon(1,:);
end

lon = mod(lon, 360);
lat = mod(lat, 360);

%% Latitude to interval [-90,90]

% first to interval [-180, 180]
ind = lat>180;
lat(ind) = lat(ind) - 360;
% then to interval [-90, 180]
ind = lat<-90;
lon(ind) = mod(lon(ind)+180, 360);
lat(ind) = -180 - lat(ind);
% and finally to interval [-90,90]
ind = lat>90;
lon(ind) = mod(lon(ind)+180, 360);
lat(ind) = 180 - lat(ind);

%% Longitude to interval [-180, 180]

lon(lon>180) = lon(lon>180) - 360;

if nargin == 1 && nargout == 1
  lon = [lon(:)'; lat(:)'];
end