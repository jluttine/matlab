
function [D, dD1, dD2] = dist_earth(coord1, coord2)
% [D, dD] = dist_earth(coord1, coord2)
% coord1 is 2 x N
% coord2 is 2 x M
% Returns N x M matrix D of mutual distances.

if nargin < 2 || isempty(coord2)
  coord2 = coord1;
end

R = 6371.01; % Spherical Earth radius approximation

n1 = cols(coord1);
n2 = cols(coord2);

D = zeros(n1, n2);

coord1 = pi/180 * coord1;
coord2 = pi/180 * coord2;
if nargout >= 2
  dD1 = zeros([n1, n2, 2]);
end
if nargout >= 3
  error('Hmm.. maybe you shouldn''t use third output? :)');
  dD2 = zeros([n1, n2, 2]);
end

for i=1:n1
  lat1 = coord1(2,i);
  lat2 = coord2(2,:);
  dlon = (coord2(1,:) - coord1(1,i));
  f = sin(lat2).*sin(lat1) + cos(lat2).*cos(lat1).*cos(dlon);
  f(f>=1-eps) = 1-eps; % correction because of numerical errors, f should be [-1,1]
  f(f<=eps-1) = eps-1; % correction because of numerical errors, f should be [-1,1]
  D(i,:) = R * acos(f);
  if ~isreal(D)
    coord1
    f
%    D
    error('Oohps! Complex distance..');
  end
  %dD_df(isinf(dD_df)) = 0;%-R ./ sqrt(1-f.^2);
  
  if nargout >= 2
    dD_df = -R ./ sqrt(1-f.^2);
    dD1(i,:,1) = pi/180 * dD_df .* cos(lat1) .* cos(lat2) .* sin(dlon);
    dD1(i,:,2) = pi/180 * dD_df .* (cos(lat1).*sin(lat2) - sin(lat1).* ...
                                    cos(lat2).*cos(dlon));
  end
% $$$   if any(isnan(dD1))
% $$$     dD1
% $$$     error('WTF?');
% $$$   end
% $$$   if nargout >= 3
% $$$     dD2(i,:,1) = -dD1(i,:,1);
% $$$     dD2(i,:,2) = dD_df .* (sin(lat1).*cos(lat2) - cos(lat1).*sin(lat2).*cos(dlon));
% $$$   end
end
