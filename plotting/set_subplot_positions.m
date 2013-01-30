
% SET_SUBPLOT_POSITIONS(HAX, ROW, COL, MAR, SEP)
%
% Places subplots to a nice layout.
%
% HAX : handles to axes
% ROW : number of rows in the layout
% COL : number of columns in the layout
% MAR : margins [left top right bottom]
% SEP : separating space between subplots [horizontal vertical]
%
% All the units are normalized.

% TODO: Maybe should take handle to figure instead of axes?
% TODO: Let the user change 'Units' to something else than 'normalized'

% Copyright (c) 2010 Jaakko Luttinen

function set_subplot_positions(hax, rows, columns, margins, separators)

% Dimensions of the area for plots
area_width = 1 - margins(1) - margins(3);
area_height = 1 - margins(2) - margins(4);

plot_width = (area_width - (columns-1) * separators(1)) / columns;
plot_height = (area_height - (rows-1) * separators(2)) / rows;

ind = 0;
for i=1:rows
  for j=1:columns
    ind = ind + 1;
    if ind <= numel(hax)
      pos(1) = margins(1) + (j-1)*(plot_width+separators(1));
      pos(3) = plot_width;
      pos(2) = 1 - (margins(2) + plot_height + (i-1)* ...
                    (plot_height+separators(2)));
      pos(4) = plot_height;
      units = get(hax(ind), 'Units');
      set(hax(ind), 'Units', 'normalized');
      set(hax(ind), 'Position', pos);
      set(hax(ind), 'Units', units);
    end
  end
end
