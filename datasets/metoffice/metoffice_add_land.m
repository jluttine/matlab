% Add land and sea-ice grid boxes to the reconstructed fields
function Yrec_out = metoffice_add_land( Yrec, dataset )

[ mask, lat, lon ] = metoffice_get_mask( dataset );

% path from system
datapath = metoffice_get_path();
load( [ datapath, dataset ], 'time' )

% Add land and ice areas to the reconstruction fields
Yrec_out = NaN( length(lat)*length(lon), length(time) );
Yrec_out( logical(mask(:)), : ) = Yrec;



function Y = metoffice_add_land_old(X, maskfile)

S = load('/share/climate/data/UK_Met_Office/RecTest/mask');
Y = nan(numel(S.mask),size(X,2));
Y(S.mask(:)==1,:) = X;
% $$$ S = load('/share/climate/data/UK_Met_Office/RecTest/land_hadgem');
% $$$ land = (S.land(:) >= 0.5);
% $$$ Y = nan(length(land),size(X,2));
% $$$ Y(~land,:) = X;
