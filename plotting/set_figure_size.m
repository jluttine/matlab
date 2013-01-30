% SET_FIGURE_SIZE - Set the size of a figure.
%
% SET_FIGURE_SIZE(WIDTH,HEIGHT,FIG,UNITS)
%
% If FIG is not given, GCF is used. The default for UNITS is 'centimeters'.

% Last modified 2010-06-11
% Copyright (c) Jaakko Luttinen (jaakko.luttinen@tkk.fi)

function set_figure_size(width,height,fig,units)

if nargin < 3
  fig = gcf;
end

if nargin < 4
  units = 'centimeters';
end

old_units = get(fig, 'Units');
set(fig, 'Units', units);
set(fig, 'PaperUnits', units);
pos = get(fig, 'Position');
set(fig, 'Position', [pos(1) pos(2) width height]);
set(fig, 'PaperPosition', [0 0 width height]);
set(fig, 'PaperSize', [width height]);
%set(fig, 'PaperPositionMode', 'auto', 'PaperSize', [width height]);
set(fig, 'Units', old_units);
