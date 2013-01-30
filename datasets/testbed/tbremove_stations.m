%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function data_struct = remove_stations(data_struct, stations)

data_struct.observations(stations,:) = [];
data_struct.stations(stations) = [];
data_struct.coordinates(stations,:) = [];
