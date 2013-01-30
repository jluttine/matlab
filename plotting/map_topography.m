% MAP_TOPOGRAPHY - Plots 5-minute resolution topography map.

% Last modified 2010-06-11
% Copyright (c) Jaakko Luttinen (jaakko.luttinen@tkk.fi)

function map_topography()

% Draw high resolution topography
values = m_etopo2('contourf', 100, 'EdgeColor', 'None');

% Set colormap
map_colormap();

% Set the range of the coloring (sea level is white).
max_value = max(values, [], 2);
set( gca, 'clim', [-max_value(1), max_value(1)]);
