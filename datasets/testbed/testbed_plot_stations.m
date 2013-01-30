% TESTBED_PLOT_STATIONS - Plots the stations in the Testbed dataset.
%
% TESTBED_PLOT_STATIONS(COORDINATES, ...)
%
% COORDINATES is a matrix that has longitudes as the first row and latitudes
% as the second row. It can also be a struct which has COORDINATES as a
% field: COORDINATES.coordinates.

% Last modified 2010-06-11
% Copyright (c) Jaakko Luttinen (jaakko.luttinen@tkk.fi)

function testbed_plot_stations(coordinates, varargin)

% Default options
options = struct(...
    'Style', 'o', ...
    'MarkerSize', 4, ...
    'MarkerEdgeColor', 'k', ...
    'MarkerFillColor', [0.7 0 0.7]);
[options, errmsg, remopts] = argparse(options, varargin{:});
error(errmsg);

if nargin < 1
  % Default set of stations
  data = testbed_loaddata();
  data = testbed_preprocess(data);
  coordinates = data.coordinates;
end

if ~isnumeric(coordinates) && isstruct(coordinates)
  coordinates = coordinates.coordinates;
end

% Plot the stations
map_plot(coordinates', ...
         options.Style, ...
         'MarkerSize',      options.MarkerSize, ...
         'MarkerEdgeColor', options.MarkerEdgeColor, ...
         'MarkerFaceColor', options.MarkerFillColor, ...
         remopts{:});

