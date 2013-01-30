function [data, folder, maskfile] = metoffice_getdata(dataset, anomalies, remove_val)

if nargin < 3 || isempty(remove_val)
    remove_val = 1;
end

%
% Load the data and remove land areas
%

logname = getenv('LOGNAME');
if strcmp( logname, 'alexilin' ) || strcmp( logname, 'hadia' )
    jaakkopath
end

os = system_dependent('getos');
infinland = ~isempty(strfind(os,'Ubuntu')) || ~isempty(strfind(os,'Linux'));
if infinland
    datapath = '/share/climate/data/UK_Met_Office/RecTest/DATASETS';
    % Folder for saving the results
    if strcmp( logname, 'alexilin' )
        folder = '/share/climate/data/UK_Met_Office/RecTest/RESULTS';
    else
        folder = '/share/work/jluttine/metoffice/rectest';
    end
else
    datapath = '/net/home/h03/hadia/data/RecTest/';
    folder = '/net/home/h03/hadia/data/RecTest/RESULTS';
end

switch dataset
case {'hadsst2d1','hadsst2d2','mohsst5d1','mohsst5d2'}
    datafile = [ datapath '/' dataset ];
    maskfile = 'mask';
case {'aari','aari8'}
    datafile = [ datapath '/' dataset ];
    maskfile = 'mask_seaice';
otherwise
    error(sprintf('Unknown dataset %s',dataset));
end

data = load(datafile);

[LON,LAT] = meshgrid(data.lon,data.lat);
data.coordinates = double(metoffice_remove_bins([LON(:),LAT(:)],maskfile))';

if isempty( findstr( maskfile, 'seaice' ) )
    
    % SST
    Itrain = ~isnan(data.sst);
    if remove_val
        Ival = get_valset( remove_val, data.time, Itrain );
        data.sst(Ival) = NaN;
        fprintf( 'Validation set %d removed\n', remove_val )
    end

    % Remove land areas
    data.data = double(metoffice_remove_bins(data.sst,maskfile));
    data = rmfield( data, 'sst' );

else

    % Sea ice
    % Remove areas
    data.data = double(metoffice_remove_bins(data.sice,maskfile));
    data = rmfield( data, 'sice');
    
end

data.clim = double(metoffice_remove_bins(data.clim,maskfile));
data.gridsize = metoffice_remove_bins(cosd(LAT(:)),maskfile);

if anomalies==true
  % Remove climatological monthly averages
  dv = datevec( data.time ); month = dv(:,2);
  
  data.data = data.data - data.clim(:,month);
  clear month dv

end
