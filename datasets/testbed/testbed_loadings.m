% TESTBED_LOADINGS - Plots interpolated spatial components.
% 
% Plots interpolated spatial components Q.W Q on map according to the given
% coordinates Q.coordinates.
%
% HAX = TESTBED_LOADINGS(Q,...)
%
% Optional parameters:
%
% 'resolution' : vector of longitudinal and latitudinal resolution
% 'components' : indeces of the components to plot
% 
% also, the optional parameters of MAP_INTERPOLATE may be used.
%
% See also MAP_INTERPOLATE.

% Last modified 2010-06-11
% Copyright (c) Jaakko Luttinen (jaakko.luttinen@tkk.fi)

function hax = testbed_loadings(W, coordinates, varargin)

options = struct( 'resolution', [60 40]);

% Check arguments
[options, errmsg, remopts] = argparse(options, varargin{:});
error(errmsg);

figure
[lon,lat] = testbed_grid(options.resolution(1), options.resolution(2));
hax = map_interpolate(coordinates(:,1:2)', ...
                      W, ...
                      lon, lat, ...
                      remopts{:});

%hcb = colorbars(hax, 'SouthOutside');
testbed_coast(hax);
%map_grid(hax, 'ytick', [], 'xtick', []);

if nargout < 1
  clear hax;
end