
% DATA = HADSST2_PARSEDATA()

% Last modified 2011-01-11
% Copyright (c) Jaakko Luttinen (jaakko.luttinen@tkk.fi)

function data = hadsst2_parsedata

filename = '/share/climate/data/UK_Met_Office/HadSST2/HadSST2_nobs.nc';

longitude = nc_varget(filename, 'longitude');
latitude = nc_varget(filename, 'latitude');
[LON,LAT] = meshgrid(longitude, latitude);
coordinates = [LON(:)';LAT(:)'];

time = nc_varget(filename, 'time');
time = time(:)' + datenum([1826 1 1]);

observations = nc_varget(filename, 'data');
%observations(observations==raw.data.FillValue) = NaN;
observations = reshape(observations, [length(time), size(coordinates, 2)])';

data = struct('observations', observations, ...
              'coordinates', coordinates, ...
              'longitude', longitude, ...
              'latitude', latitude, ...
              'time', time);
