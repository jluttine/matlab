% Loads the true SST fields from the HadGem1 runs, the input
% is the name of the file name with the training data
%
function [ true_sst, Itest, Itrain, time, train_sst ] = ...
    load_testdata( dataset )

switch dataset
case 'hadsst2d1'
    datafile = 'run1_ts_A1_5x5';
case 'hadsst2d2'
    datafile = 'run2_ts_A1_5x5';
end

% Load training data
datapath = metoffice_get_path();
train = load( [ datapath, dataset ], 'time', 'sst' );
%train = load( [ datapath '/RecTest/DATASETS/' dataset ], 'time', 'sst' );
Itrain = ~isnan( train.sst );

% Load test data
testpath = metoffice_get_testpath();
hadgem = load( [ testpath, datafile ], 'sst', 'time' );
% $$$ hadgem = load( [ datapath '/HadGem1/' datafile ], 'sst', 'time' );

[ time, Ihadsst, Ihadgem ] = intersect( train.time, hadgem.time );

% Take the time period included in the test data
true_sst = hadgem.sst(:,Ihadgem);

% Load mask
[ mask, lat, lon ] = metoffice_get_mask( dataset );

Itest = logical( ones(size(true_sst)) );
Itest( ~mask(:), : ) = 0;
Itest( Itrain ) = 0;
% Possible missing data in model runs
Itest( isnan(true_sst) ) = 0;
% Remove ice from test data
Iice = true_sst < -1.8;
Itest(Iice) = 0;

% Remove Caspian sea (detached area with no observations)
mask_Caspian = zeros(length(lat),length(lon));
%          lat  lon
ilat = find(lat==47.5);
ilon = find(lon==52.5);
mask_Caspian( ilat, ilon ) = 1;
mask_Caspian( ilat+1, ilon ) = 1;
mask_Caspian( ilat+2, ilon ) = 1;
mask_Caspian( ilat, ilon-1 ) = 1;
mask_Caspian = reshape( mask_Caspian, length(lat)*length(lon), 1 );
%lmap( mask_Caspian, lat, lon );

Itest( logical(mask_Caspian), : ) = 0;

%true_sst(~(Itest | Itrain)) = NaN;

if nargout >= 5
    train_sst = train.sst(:,Ihadsst);
end
