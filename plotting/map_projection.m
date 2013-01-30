% MAP_PROJECTION - Set the projection
%
% 'GLOBAL-ELLIPSE'
% 'GLOBAL-RECT'
% 'TESTBED'

function map_projection(projection)

if nargin < 1
  projection = 'global-ellipse';
end

switch lower(projection)
 case 'global-ellipse'
  m_proj('mollweide', 'clongitude', 0); % straight latitudes (equal area)
  %m_proj('hammer-aitoff', 'clon', 0); % curved latitudes (equal area)
 case 'global-rect'
  m_proj('miller');
 case 'testbed'
  m_proj('Mercator', 'lon',[22.5,26.8], 'lat',[59.74,61.0])
 otherwise
  error('Unknown projection')
end

