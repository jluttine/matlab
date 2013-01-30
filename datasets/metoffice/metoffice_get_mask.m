% Get the mask defining the points where the reconstructions should be
% computed.
function [ mask, lat, lon ] = metoffice_get_mask( dataset )

% filename from dataset
switch dataset
case {'hadsst2d1','hadsst2d2','mohsst5d1','mohsst5d2'}
    maskfile = 'mask';
case {'aari','aari8'}
    maskfile = 'mask_seaice';
otherwise
    error(sprintf('Unknown dataset %s',dataset));
end

% path from system
datapath = metoffice_get_path();

% load mask
if nargout == 3
    load( [ datapath, maskfile ], 'mask', 'lat', 'lon' )
elseif nargout == 1
    load( [ datapath, maskfile ], 'mask' )
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function m = metoffice_get_mask_old(maskfile)

if nargin == 0
  warning('Give the maskfile');
  maskfile = 'mask';
end

os = system_dependent('getos');

if ~isempty(strfind(os,'Ubuntu'))
  % Finland
  maskfile = [ '/share/climate/data/UK_Met_Office/RecTest/DATASETS/' maskfile ];
elseif ~isempty(strfind(os,'Linux'))
  % ICS GridEngine
  maskfile = [ '/share/climate/data/UK_Met_Office/RecTest/DATASETS/' maskfile ];
else
  % MetOffice
  maskfile = [ '/net/home/h03/hadia/data/RecTest/' maskfile ];
end

S = load(maskfile);
m = S.mask(:)~=0;
