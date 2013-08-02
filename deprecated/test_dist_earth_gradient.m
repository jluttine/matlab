function test_dist_earth_gradient

%rand('state', 1);

n = 100;
lon = 5 * rand(1,n);
lat = 5 * rand(1,n);
% $$$ lon = 360 * rand(1,n);
% $$$ lat = 180 * rand(1,n);

coord = [lon; lat];

coordT = coord';

mycheckgrad(@func, coordT(:), 1e-2);

%coordTmin = coordT;
coordTmin = minimize(coordT(:), @func, 150);
coordmin = reshape(coordTmin, [n,2])';

figure
mapproj('global-ellipse');

[lon, lat] = lonlat2standard(coordmin(1,:), coordmin(2,:));
mapplot(lon, lat, 'r+');
hold on
mapcoast
%error('j')


function [f, df] = func(x)

n = length(x) / 2;
coord = reshape(x, [n,2])';

%[D, dD1] = dist_earth(coord, [0:10:50; 10:20:110]);

% Square root helps in placing the locations uniformly
[D, dD1] = dist_earth(coord, coord);
f = -sum(sqrt(D(:))) / n;

%df = colsum(D,dD1);
%df = -2 * df(:) / n;
df = 0.5*colsum(bsxfun(@times,1./sqrt(D),dD1));
df = -2 * df(:) / n;

% $$$ f = D(1,1)
% $$$ df = dD1(1,1,:);
% $$$ df = [df(:); zeros((n-1)*2,1)];

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5

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