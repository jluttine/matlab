function data = tbsubset(data, stations, times)

remove = setdiff(1:rows(data.observations), stations);
data = tbremove_stations(data, remove);
data.time = data.time(times);
data.observations = data.observations(:,times);