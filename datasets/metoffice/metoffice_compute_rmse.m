% Usage: compute_rmse( Yrec, dataset, remove_val )
%        compute_rmse( name_of_file_with_reconstructions )
function [ out, lat, lon, time ] = metoffice_compute_rmse( Yrec, dataset, remove_val )

if isstr(Yrec)
    setup = get_setup( Yrec )

    % Load reconstructions
    datapath = metoffice_get_path();
    load( [ datapath, Yrec ], 'Yrec' );
    %load( [ datapath '/RecTest/RESULTS/' Yrec ], 'Yrec' );
else
    setup.dataset = dataset;
    setup.remove_val = remove_val;
    setup
end

% Remove values below -1.8
Yrec( Yrec < -1.8 ) = -1.8;

[ true_sst, Itest, Itrain, time, train_sst ] = metoffice_get_testdata( setup.dataset );
%[ true_sst, Itest, Itrain, time, train_sst ] = load_testdata( setup.dataset );

if setup.remove_val
    Ival = get_valset( setup.remove_val, time, Itrain );
    Itrain = Itrain & ~Ival;
end

% Load lat and lon to compute the area weights
datapath = metoffice_get_path();
load( [ datapath, dataset ], 'lat', 'lon' )
%load( [ datapath '/RecTest/DATASETS/' dataset ], 'lat', 'lon' )
lts = repmat( lat, 1, length(lon) );
weights = cosd( lts(:) );
clear lts lat lon

err = true_sst - Yrec;
out.train = get_rmse( err, Itrain, weights );
out.test = get_rmse( err, Itest, weights );
if setup.remove_val
    out.val = get_rmse( err, Ival, weights );
end

%out.Yrec = Yrec;

out.trdata = get_rmse( train_sst-Yrec, Itrain, weights );

%out.test.err = err;
%out.test.err(~Itest) = NaN;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = get_rmse( err, I, weights )

err2 = err.^2;
err2 = bsxfun( @times, weights, err2 );
Iw = bsxfun( @times, weights, I );

err2(~I) = 0;
out.rmse_fld = sqrt( sum(err2,2)./sum(I,2) );

out.rmse_ts = sqrt( sum(err2,1)./sum(I,1) );

out.rmse = sqrt( sum( err2(:) )/sum(I(:)) );

if 1
    out.rmse2 = sum(err2,1)./sum(Iw,1);
    out.rmse2( isnan(out.rmse2) ) = [];
    out.rmse2 = sqrt( mean( out.rmse2 ) );
end

out.err_fld = sqrt(sum(err2,2));
out.err_ts = sqrt(sum(err2,1));

% Smoothed time series: take 5 year moving average
err2_ts = sum(err2,1);
nobs_ts = sum(I,1);

nfilt = 5*12;
T = size(err2,2);
err2_sm = conv( ones(1,nfilt), err2_ts );
err2_sm = err2_sm( nfilt/2+(1:T) );
nobs_sm = conv( ones(1,nfilt), nobs_ts );
nobs_sm = nobs_sm( nfilt/2+(1:T) );

out.rmse_smts = sqrt( err2_sm./nobs_sm );

out.err_smts = sqrt(err2_sm);

% Median
out.medwe = median( err(I) );
out.medwe_fld = NaN*zeros(size(err,1),1);
for i = 1:size(err,1)
    out.medwe_fld(i) = median( err(i, I(i,:) ) );
end

out.medwe_ts = NaN*zeros(1,size(err,2));
for i = 1:size(err,2)
    out.medwe_ts(i) = median( err( I(:,i), i ) );
end
