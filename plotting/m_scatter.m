function han = m_scatter(lon, lat, varargin)
% function m_scatter(lon, lat, S, C)

m_proj('miller','lat',82);
%m_grid('linestyle','none','box','on','tickdir','out');
[x,y] = m_ll2xy(lon, lat);
%drawnow
h = scatter(x, y, varargin{:});
line(x(1:10), y(1:10))
%plot(x,y,'.')
%get(gca, 
z = findobj(gca, 'type', 'scattergroup')
m_grid('linestyle','none','box','on','tickdir','out')
%h = scatter(x, y, varargin{:});

if nargout == 1
  han = h;
end
