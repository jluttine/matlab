% Get validation set indices for SST tests.
function Ival = get_valset( remove_val, time, Itrain )

if nargin == 1
    datafile = 'test_hadsst2_d1.mat';

    warning( sprintf( 'Load validation set for %s', datafile ) )
    
    load( [ datapath '/RecTest/' datafile ],...
          'sst', 'time', 'lat', 'lon', 'clim' )
    
    Itrain = ~isnan( sst );
end

switch remove_val
case 1
    % The validation set is created from the recent data using the
    % flipped pattern from the first 100 years with the interval of 13
    % months.
    
    I = 1:13:(12*100);
    Tval = length(I);

    % Take the last period with the selected coverage as validation test
    
    Ival = logical(zeros(size(Itrain)));
    Ival(:,1:Tval) = ~Itrain(:,I);
    
    % Revert (flip) the order
    Ival = fliplr(Ival);
    Ival = Ival & Itrain;

case 2
    % The same as 1 but the interval is 11 months

    I = 12:11:(12*100);
    Tval = length(I);
    
    % Take the last period with the selected coverage as validation test

    Ival = logical(zeros(size(Itrain)));
    Ival(:,1:Tval) = ~Itrain(:,I);

    % Revert (flip) the order
    Ival = fliplr(Ival);
    Ival = Ival & Itrain;
    
case 3
    % Take the flipped coverage from the past

    Ival = fliplr( ~Itrain );
    Ival = Ival & Itrain;
    
case 4
    
    % Take the flipped coverage from the first 100 years
    
    I100 = Itrain(:,1:12*100);
    I100 = fliplr(I100);
    Ival = logical(zeros(size(Itrain)));
    Ival(:,1:12*100) = I100 & Itrain(:,1:12*100);
    
case {11,12,13,14,15,16,17,18,19,20}
    % The validation set is created using a coverage from a 10-year
    % period from the past. The value of remove_val (from 11 to 20)
    % defines the 10-year period used for creating the validation set.
    
    ix = remove_val - 10;
    % Coverage interval for the validation data
    tstart = 1860:10:1950;
    tstart = datenum( [ tstart(ix) 01 15 ] );
    ix = find( time == tstart );
    Icovval = ix:ix+12*10-1; % 10 year interval
    
    Ival = logical(zeros(size(Itrain)));
    Ival(:,1:size(Icovval,2)) = ~Itrain(:,Icovval);
    % Revert the order
    Ival = fliplr(Ival);
    Ival = Ival & Itrain;

otherwise
    error( 'Unknown remove_val' )
    
end

if nargout == 0
    
    sst_train2 = sst;
    sst_train2(Ival) = NaN;
    
    Itrain2 = Itrain & ~Ival;
    
    plot_ts( sum(Itrain2,1), 'time', time )
    plot_ts( sum(Ival,1), 'time', time )
    
    mapts( sst_train2, lat, lon, time )

    sst_train2 = sst;
    sst_train2(~Ival) = NaN;
    mapts( sst_train2, lat, lon, time )

    fprintf( 'Number of training points: %d\n', sum(Itrain(:)) )
    fprintf( 'Number of validation points: %d\n', sum(Ival(:)) )
    
    clear Ival
end