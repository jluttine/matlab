% COLORBAR_POSITION - Set the position of a colorbar.
%
% HCB = COLORBAR_POSITION(SIZE,SEP,LOCATION,HAX,...)
%
% SIZE     : desired size of the colorbar
% SEP      : separation between the axes and the colorbar
% LOCATION : location relative to the axes [default: 'EastOutside'] (see COLORBAR)
% HAX      : handle to parent axes [default: GCA]
%
% For optional parameters, see COLORBAR.
%
% HCB : handle to the colorbar object

% Last modified 2010-06-11
% Copyright (c) Jaakko Luttinen (jaakko.luttinen@tkk.fi)

function hcb = colorbar_position(sz,sep,location,hax,varargin)

error('Deprecated. Use COLORBARS.');

if nargin < 3
  location = 'EastOutside';
end

if nargin < 4
  hax = gca;
end


hcb = zeros(size(hax));

for d=1:numel(hax)
  
  % Get axes position
  pos_ax = get(hax(d), 'Position');
  
  % Create colorbar object
  hcb(d) = colorbar('peer', hax(d), location, varargin{:});
  
  % Use the same units
  units = get(hax(d), 'Units');
  set(hcb, 'Units', units);
  
  pos = get(hcb, 'Position');

  if strncmpi(location,'North',5) || strncmpi(location,'South',5)
    pos(4) = sz;
  else
    pos(3) = sz;
  end

  if strcmpi(location,'North')
    pos(2) = pos_ax(2) + pos_ax(4) - pos(4) - sep;
  elseif strcmpi(location,'South')
    pos(2) = pos_ax(2) + sep;
  elseif strcmpi(location,'East')
    pos(1) = pos_ax(1) + pos_ax(3) - pos(3) - sep;
  elseif strcmpi(location,'West')
    pos(1) = pos_ax(1) + sep;
  elseif strcmpi(location,'NorthOutside')
    pos(2) = pos_ax(2) + pos_ax(4) + sep;
  elseif strcmpi(location,'SouthOutside')
    pos(2) = pos_ax(2) - sep - pos(4);
  elseif strcmpi(location,'EastOutside')
    pos(1) = pos_ax(1) + pos_ax(3) + sep;
  elseif strcmpi(location,'WestOutside')
    pos(1) = pos_ax(1) - sep - pos(3);
  else
    error('Unknown location given')
  end
  
  set(hcb(d), 'Units', units, 'Position', pos);

end

if nargout == 0
  clear hcb;
end