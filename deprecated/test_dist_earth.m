function test_dist_earth
warning('This function is deprecated')


% Grid
% $$$ lon = linspace(22,28,100);
% $$$ lat = linspace(55,65,60);
lon = linspace(-180,180,100);
lat = linspace(-90,90,60);
[LAT,LON] = meshgrid(lat,lon);

% Distances from..
coord0 = [25; 60];
D = dist_earth([LON(:), LAT(:)]', coord0);
D = reshape(D, [length(lon), length(lat)]);

% Plot
figure
mapproj('global-ellipse')
% $$$ mapproj('testbed')
mapcolor(lon, lat, D);
mapcoast;