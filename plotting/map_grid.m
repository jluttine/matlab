
% MAP_GRID
%
% Draws the map grid and box using M_MAP package.
%
% See M_COAST and M_USERCOAST for optional arguments.
%
% By default, 'box' is set to 'on'. If 'hammer-aitoff' projection is
% used, 'xtick' is set to [] by default.

% Last modified 2010-12-17
% Copyright (c) Jaakko Luttinen (jaakko.luttinen@tkk.fi)

function map_grid(varargin)

global MAP_PROJECTION

%opts = struct('xtick', [], 'xticklabels', []);

% Some default parameters
opts.box = 'on';
if strcmpi(MAP_PROJECTION.name, 'hammer-aitoff')
  opts.xtick = [];
  %opts.xtick = 'xtick', [];
end

% Read user's parameters
if nargin >= 1 && isnumeric(varargin{1})
  hax = varargin{1};
  [opts, errmsg, remopts] = argparse(opts, varargin{2:end});
else
  hax = NaN;
  [opts, errmsg, remopts] = argparse(opts, varargin{:});
end
opts = struct2list(opts);
params = {opts{:}, remopts{:}};

% Draw grid
if ~isnan(hax)
  old_axes = gca();
  for n=1:numel(hax)
    axes(hax(n));
    m_grid(params{:});
  end
  axes(old_axes);
else
  m_grid(params{:});
end
