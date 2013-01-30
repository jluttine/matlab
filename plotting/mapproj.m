
function mapproj(projection)

if nargin < 1
  projection = 'global';
end

switch lower(projection)
 case 'global-ellipse'
  m_proj('hammer-aitoff', 'clon', 0);
 case 'global-rect'
  m_proj('miller');
 case 'testbed'
  m_proj('Mercator', 'lon',[22.5,26.8], 'lat',[59.74,61.0])
end

