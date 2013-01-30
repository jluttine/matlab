% MAP_INTERPOLATE - Plots interpolated values over a map.
%
% Plots a smooth coloring over a map based on interpolation using values at
% some specific locations. It is possible to make several interpolated plots
% if the set of locations is identical.
%
% MAP_INTERPOLATE(C, Y, LONH, LATH, ...)
%
% where 
%
% C    : 2 x N matrix of longitudes and latitudes
% Y    : N x K matrix of values
% LONH : vector defining the longitudinal grid points
% LATH : vector defining the latitudinal grid points
%
% Optional parameters:
%
% 'method'   : 'nn' is nearest neighbour interpolation
%              'rbf' is linear radial basis functions interpolation (default)
% 'layout'   : vector of the number of rows and columns of the subplot
%              layout (default is square layout)
% 'distfunc' : distance function to be used (default: 'haversine')

% Last modified 2010-06-11
% Copyright (c) Jaakko Luttinen (jaakko.luttinen@tkk.fi)

function hax = map_interpolate(C, Y, lonh, lath, varargin)

options = struct( ...
    'method',   'rbf', ...
    'layout',   [],   ...
    'distfunc', 'haversine');

% Check arguments
[options, errmsg] = argparse(options, varargin{:});
error(errmsg)

% Form a regular grid
[LONH,LATH] = meshgrid(lonh,lath);

% Evaluate distances for interpolation
Ch = [LONH(:),LATH(:)]';
D = dist_haversine(C, C);
Dh = dist_haversine(Ch, C);

% Settings for multiple plots
K = size(Y,2);
if K > 1 && ~isempty(options.layout)
  M = options.layout(1);
  N = options.layout(2);
else
  M = ceil(sqrt(K));
  N = ceil(K/M);
end
hax = zeros(K,1);

% Plot maps
for k=1:K
  
  hax(k) = subplot(M,N,k);

  % Interpolate
  yh = interpolate(D, Y(:,k), Dh, 'method', options.method);
  yh = reshape(yh, size(LONH));

  % Plot
  map_pcolor(LONH, LATH, yh);

end

if nargout < 1
  clear hax;
end
