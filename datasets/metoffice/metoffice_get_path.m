function datapath = metoffice_get_path()

os = system_dependent('getos');
if ~isempty(strfind(os,'Ubuntu'))
  % Finland
  datapath = '/share/climate/data/UK_Met_Office/RecTest/DATASETS/';
elseif ~isempty(strfind(os,'Linux'))
  % ICS GridEngine
  datapath = '/share/climate/data/UK_Met_Office/RecTest/DATASETS/';
else
  % MetOffice
  datapath = '/net/home/h03/hadia/data/RecTest/';
end
