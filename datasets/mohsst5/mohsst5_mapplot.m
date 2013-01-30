% mohsst5_mapplot(w, data)

function hax = mohsst5_mapplot(w, varargin)

options = struct('colormap', 'centered', ...
                 'data', [], ...
                 'ytick', []);

[options, errmsg] = argparse(options, varargin{:});
error(errmsg);

if nargin < 2 || isempty(options.data)
  data = mohsst5_loaddata();
else
  data = options.data;
end

if ~isvector(w)
  error('The loadings w should be a vector');
end

M = size(data.observations,1);

lands = colsum(~isnan(data.observations)) == 0;
if length(w) == 1727
  % Fill the land areas to the matrix
  tmp = nan(M,1);
  tmp(~lands) = w;
  w = tmp;
end

map_projection('global-ellipse');
hax = map_pcolor(data.longitude, ...
                 data.latitude, ...
                 reshape(w, [length(data.longitude) length(data.latitude)]));

switch options.colormap
 case 'centered'
  colormap_redblue();
  clim_centered();
 case 'positive'
  colormap_scale();
  clim_positive();
end

map_coast();
map_grid('yticklabels', [], 'xtick', [], 'ytick', options.ytick);

if nargout < 1
  clear hax;
end

