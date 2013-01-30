function ind = mohsst5_points_to_grid_index(C)

%data = mohsst5_loaddata();

lon = C(1,:);
lat = C(2,:);

N_lon = 72;
N_lat = 36;

% Random locations to 5x5-grid
ind_lon = floor(lon/5)+(N_lon/2+1);
ind_lat = N_lat - (floor(lat/5)+(N_lat/2+1)) + 1;

ind = (ind_lon-1)*N_lat + ind_lat;