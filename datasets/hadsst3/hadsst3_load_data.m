
% DATA = HADSST2_PARSEDATA()

% Last modified 2011-01-11
% Copyright (c) Jaakko Luttinen (jaakko.luttinen@tkk.fi)

function data = hadsst3_load_data()

filename = '/share/climate/data/HadSST3/HadSST3_mostlikely.nc';
%filename = '/share/climate/data/UK_Met_Office/HadSST2/HadSST2_nobs.nc';

longitude = nc_varget(filename, 'longitude');
latitude = nc_varget(filename, 'latitude');
[LON,LAT] = meshgrid(longitude, latitude);
coordinates = [LON(:)';LAT(:)'];

time = nc_varget(filename, 'time');
time = time(:)' + datenum([1850 0 0]);

observations = nc_varget(filename, 'sst');
%observations(observations==raw.data.FillValue) = NaN;
observations = reshape(observations, [length(time), size(coordinates, 2)])';

data = struct('observations', double(observations), ...
              'coordinates', double(coordinates), ...
              'longitude', double(longitude), ...
              'latitude', double(latitude), ...
              'time', double(time));
