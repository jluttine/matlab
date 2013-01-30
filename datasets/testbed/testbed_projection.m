% TESTBED_PROJECTION - Set the map projection for Testbed plots.
%
% [MINLON,MAXLON,MINLAT,MAXLAT] = TESTBED_PROJECTION()

% Last modified 2010-06-11
% Copyright (c) Jaakko Luttinen (jaakko.luttinen@tkk.fi)

function [minlon, maxlon, minlat, maxlat] = testbed_projection()

minlon = 22.5;
maxlon = 26.8;
minlat = 59.75;
maxlat = 61.0;
m_proj('Mercator', 'lon',[minlon,maxlon], 'lat',[minlat,maxlat]);

if nargout == 0
  clear minlon maxlon minlat maxlat;
end
