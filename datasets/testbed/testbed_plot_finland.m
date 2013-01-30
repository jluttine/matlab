% TESTBED_PLOT_FINLAND - Plot the map of Finland with a red rectangle
% showing the Testbed research area.

% Last modified 2010-06-11
% Copyright (c) Jaakko Luttinen (jaakko.luttinen@tkk.fi)

function testbed_plot_finland()

% Get the research area
[minlon,maxlon,minlat,maxlat] = testbed_projection();

% Set the projection
m_proj('Mercator', 'lon',[19,32], 'lat',[59.3,70.2])

% Plot the country boundary
m_plotbndry('finland', 'color', 'k');

% Plot the research area
m_line([minlon,minlon,maxlon,maxlon,minlon], ...
       [minlat,maxlat,maxlat,minlat,minlat], ...
       'color', 'r', 'LineWidth', 2);
m_grid('box','off', 'linestyle','none', 'xticklabels',[], 'yticklabels',[], ...
       'xtick',[], 'ytick',[]);
