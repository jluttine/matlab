% POINTS_ON_SPHERE - Approximately uniformly spaced points on a sphere.
%
% Usage:
%
%   [LON, LAT] = POINTS_ON_SPHERE(N)
%
% N is the number of uniformly spaced latitudes. Each of these latitudes is
% divided into uniformly spaced longitudinal points such that the distance
% between the longitudinal points is approximately the same as between
% the latitudes.

% where D is the approximate distance between latitudes and points on the
% same latitude. It is also approximately the maximum distance from any
% location to the nearest point.

function [LON, LAT] = points_on_sphere(N_lat)

% Equally spaced latitudes:
% N=1: 1/2
% N=2: 1/4 3/4
% N=3: 1/6 3/6 5/6
% N=4: 1/8 3/8 5/8 7/8
% ...

% N_lat = ceil(1+pi/d);
lat = 180 * ((1:N_lat)*2-1)/(2*N_lat) - 90;
d = pi/N_lat;

% Compute the number of longitude points at each latitude
N_lon = zeros(N_lat,1);
%lon = cell(N_lat,1);
for n=1:N_lat
  % Radius of the latitude circle
  r = cosd(lat(n));
  % Circumference
  c = 2*pi*r;
  % Number of points
  N_lon(n) = ceil(c/d);
  %
  %lon{n} = 
end

N = sum(N_lon);
LON = zeros(N,1);
LAT = zeros(N,1);

n_tot = 0;
for n=1:N_lat
  ind = n_tot + (1:N_lon(n));
  LAT(ind) = lat(n);
  LON(ind) = 360 * ((1:N_lon(n))*2-1)/(2*N_lon(n)) - 180;
  n_tot = n_tot + N_lon(n);
end

if nargout < 2
  LON = [LON(:), LAT(:)]';
end