% MOHSST5_WEIGHTS - Returns weights for the spatial locations in MOHSST5
%                   dataset.
%
% Polar regions have zero weight and equator regions unit weight.
%
% The weight is the cosine of the latitude. It can be used to scale the
% noise precision or divide the noise variance.

function W = mohsst5_weights()

data = mohsst5_loaddata();

[LON,LAT] = meshgrid(data.longitude,data.latitude);

W = cosd(LAT(:));

