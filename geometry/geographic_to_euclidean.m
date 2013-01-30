% GEOGRAPHIC_TO_EUCLIDEAN - Transforms geographic coordinates to Euclidean
%
%   X = GEOGRAPHIC_TO_EUCLIDEAN(L, R)
%
% L is a 2xN matrix of geographic coordinates (longitude,latitude).
% R is the radius (default: spherical Earth radius 6371 kilometers)
% X is a 3xN matrix of Euclidean coordinates (x,y,z)

% Last modified 2010-12-17
% Copyright (c) Jaakko Luttinen (jaakko.luttinen@tkk.fi)

function X = geographic_to_euclidean(L, R)

if nargin < 2
  % Earth radius
  R = 6371;
end

X = zeros(3, size(L,2));

X(1,:) = R .* cosd(L(2,:)) .* cosd(L(1,:));
X(2,:) = R .* cosd(L(2,:)) .* sind(L(1,:));
X(3,:) = R .* sind(L(2,:));

