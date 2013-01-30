% MAP_PLOT - Plots on a map.
%
% [...] = MAP_PLOT(C, ...)
% 
% where C(1,:) is the longitudes, C(2,:) is the latitudes. For more
% information, see M_PLOT.

% Last modified 2010-06-11
% Copyright (c) Jaakko Luttinen (jaakko.luttinen@tkk.fi)

function varargout = map_plot(C, varargin)

varargout{:} = m_plot(C(1,:), C(2,:), varargin{:});
