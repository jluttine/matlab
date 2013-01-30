% MAP_PCOLOR - Plots a spatial function colored on a grid
%
% MAP_PCOLOR(LON, LAT, Z, ...)

function hax = map_pcolor(lon, lat, Z, varargin)

if isvector(lon) && isvector(lat) && numel(Z)==length(lon)*length(lat)
  [lon, lat] = meshgrid(lon,lat);
end

hax = m_pcolor(lon,lat,reshape(Z,size(lon)),varargin{:});
shading('flat');
%clim_centered();
%map_colormap();

if nargout < 1
  clear hax;
end

%shading('interp');
