function e = mohsst5_rmsew(Y, Yh)

Z = Y-Yh;
I = ~isnan(Z);

% Weight proportional to area size
load('mohsst5_data.mat', 'latitude', 'longitude');
%[lat,lon] = meshgrid(latitude, longitude);
[lon,lat] = meshgrid(longitude, latitude);
W = repmat(cosd(lat(:)), [1, size(Y,2)]);

e = rmsew(Z(I), W(I));
