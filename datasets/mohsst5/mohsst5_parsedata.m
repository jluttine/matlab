
% [DATA, DATA_RAW] = MOHSST5_PARSEDATA()

% Last modified 2011-01-11
% Copyright (c) Jaakko Luttinen (jaakko.luttinen@tkk.fi)

function [data, data_raw] = mohsst5_parsedata

filename = ['/share/climate/data/GOSTA_atlas7_MOHSST5_ssta/' ...
            'MOHSST5_ssta.cdf' ];
data = nc_getall( filename );
data_raw = data;

x = double(data.ssta.data);
x = shiftdim(x,1);
x = reshape( x, [ size(x,1)*size(x,2) size(x,3) ]);
    
tmp = sum( sum( isnan(x) ) ) / prod(size(x)) * 100;
fprintf( '%.2f%% of values are missing\n', tmp )
    
lat = double(data.Y.data);
lon = double(data.X.data);
m60 = double(data.T.data);

year = 1960 + floor(m60/12);
month = ceil(mod(m60,12));
timevec = [ year month repmat(15,length(m60),1) ];
%timevec
time = datenum( timevec )';
    
%save ssta x lat lon time
[LON, LAT] = meshgrid(lon,lat);
%[LAT, LON] = meshgrid(lat,lon);
data = struct('observations',x, 'coordinates',[LON(:)';LAT(:)'], 'longitude', ...
              lon, 'latitude', lat, 'time',time);