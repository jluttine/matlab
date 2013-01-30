% MAP_COAST - Draws the coast line using M_MAP package.
%
% MAP_COAST(...)
% MAP_COAST(HAX,...)
%
% HAX : vector of axes handles (default: GCA)
%
% Optional arguments:
%
% 'color'     : color of the coast line (default: [0.3 0.3 0.3])
% 'usercoast' : custom file for loading the coast line data (default: [])
%
% See M_COAST and M_USERCOAST for more optional arguments.

% Last modified 2010-06-11
% Copyright (c) Jaakko Luttinen (jaakko.luttinen@tkk.fi)

function map_coast(varargin)

options = struct( ...
    'color',     [0.3 0.3 0.3], ...
    'usercoast', []);

% Read user's parameters
if nargin >= 1 && isnumeric(varargin{1})
  hax = varargin{1};
  [options, errmsg, remopts] = argparse(options, varargin{2:end});
  error(errmsg);
else
  hax = gca();
  [options, errmsg, remopts] = argparse(options, varargin{:});
  error(errmsg);
end

% Draw coast
old_axes = gca();
for n=1:numel(hax)
  axes(hax(n));
  if ~isempty(options.usercoast)
    m_usercoast(options.usercoast, 'color', options.color, remopts{:});
  else
    m_coast('color', options.color, remopts{:});
  end
end
axes(old_axes);
  
