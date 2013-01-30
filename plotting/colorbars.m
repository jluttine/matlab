% COLORBARS - Puts colorbars on possibly multiple plots.
%
% COLORBARS()
%
% works identically to COLORBAR().
%
% COLORBARS(HAX, ...)
%
% runs COLORBAR(...) for each of the axes in HAX.
%
% HCB = COLORBARS(...)
% HCB = COLORBARS(HAX, ...)
%
% returns the handles to the colorbar objects.
%
% Optional arguments:
%
% 'Size'       : size of the colorbar (height/width for horizontal/vertical bar)
% 'Separation' : distance between the bar and the axes
% 'Location'   : {'North','South',...} (default: 'EastOutside', see
%                COLORBAR)
%
% See also COLORBAR.

% Last modified 2010-06-11
% Copyright (c) Jaakko Luttinen (jaakko.luttinen@tkk.fi)

function hcb = colorbars(varargin)

options = struct(...
    'Size',       [], ...
    'Separation', [], ...
    'Location',   'EastOutside');
    

if nargin >= 1 && isnumeric(varargin{1})
  hax = varargin{1};
  [options, errmsg, remopts] = argparse(options, varargin{2:end});
  error(errmsg);
else
  hax = gca;
  [options, errmsg, remopts] = argparse(options, varargin{:});
  error(errmsg);
end

hcb = zeros(size(hax));
for n=1:numel(hax)

  % Get axes position
  pos_ax = get(hax(n), 'Position');
  
  % Add a colorbar
  hcb(n) = colorbar('peer', hax(n), options.Location, remopts{:});
  
  % Use the same units
  units = get(hax(n), 'Units');
  set(hcb(n), 'Units', units);
  
  %pos = get(hcb(n), 'Position');
  pos = pos_ax;

  if ~isempty(options.Size)
    if strncmpi(options.Location,'North',5) || ...
          strncmpi(options.Location,'South',5)
      pos(3) = pos_ax(3);
      pos(4) = options.Size;
    else
      pos(4) = pos_ax(4);
      pos(3) = options.Size;
    end
    set(hcb(n), 'Units', units, 'Position', pos);
  end

  if ~isempty(options.Separation)
    if strcmpi(options.Location,'North')
      pos(2) = pos_ax(2) + pos_ax(4) - pos(4) - options.Separation;
    elseif strcmpi(options.Location,'South')
      pos(2) = pos_ax(2) + options.Separation;
    elseif strcmpi(options.Location,'East')
      pos(1) = pos_ax(1) + pos_ax(3) - pos(3) - options.Separation;
    elseif strcmpi(options.Location,'West')
      pos(1) = pos_ax(1) + options.Separation;
    elseif strcmpi(options.Location,'NorthOutside')
      pos(2) = pos_ax(2) + pos_ax(4) + options.Separation;
    elseif strcmpi(options.Location,'SouthOutside')
      pos(2) = pos_ax(2) - options.Separation - pos(4);
    elseif strcmpi(options.Location,'EastOutside')
      pos(1) = pos_ax(1) + pos_ax(3) + options.Separation;
    elseif strcmpi(options.Location,'WestOutside')
      pos(1) = pos_ax(1) - options.Separation - pos(3);
    else
      error('Unknown location given')
    end
    set(hcb(n), 'Units', units, 'Position', pos);
  end
end

if nargout < 1
  clear hcb;
end