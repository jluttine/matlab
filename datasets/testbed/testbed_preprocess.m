function data = testbed_preprocess(data)

% Remove stations with no observations
I = sum(~isnan(data.observations),2) == 0;
data = testbed_remove_stations(data, I);

% Remove station levels
data = remove_station_levels(data);

% Remove outregion stations
G_MINLON = 22.5;
G_MAXLON = 26.8;
G_MINLAT = 59.7;
G_MAXLAT = 61.0;
n = size(data.coordinates, 1);
I = [];
for i=1:n
  lon = data.coordinates(i,1);
  lat = data.coordinates(i,2);
  if lon < G_MINLON || lon > G_MAXLON || lat < G_MINLAT || lat > G_MAXLAT
    I = [I; i];
  end
end
data = testbed_remove_stations(data, I);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [coordinates, I] = remove_station_levels(coordinates)
% I is the indeces of the remaining stations

if isstruct(coordinates)
  coords = coordinates.coordinates;
else
  coords = coordinates;
end

lon = coords(:, 1);
lat = coords(:, 2);
n = size(coords, 1);

I = [];

for i=1:n
  ind = find(coords(I,1)==lon(i) & coords(I,2)==lat(i));
  if isempty(ind)
    I = [I; i];
  else
    s = [I(ind), i]; %indeces with same location
    [temp, min_i] = min(coords(s,3));
    I(ind) = s(min_i);
  end
end

if isstruct(coordinates)
  I_rm = 1:n;
  I_rm(I) = [];
  coordinates = testbed_remove_stations(coordinates, I_rm);
else
  coordinates = coordinates(I, :);
end
