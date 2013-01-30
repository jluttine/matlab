function datapath = metoffice_get_testpath()

os = system_dependent('getos');
if ~isempty(strfind(os,'Ubuntu'))
  % Finland
  datapath = '/share/climate/data/UK_Met_Office/HadGem1/';
elseif ~isempty(strfind(os,'Linux'))
  % ICS GridEngine
  datapath = '/share/climate/data/UK_Met_Office/HadGem1/';
else
  % MetOffice
  datapath = '/net/home/h03/hadia/data/RecTest/';
end
