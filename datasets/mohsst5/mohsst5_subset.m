function data = mohsst5_subset(data, stations, times)

[M,N] = size(data.observations);

if numel(stations) == 1 && stations < 1
  permM = randperm(M);
  subM = ceil(stations * M);
  stations = permM(1:subM);
  stations = sort(stations);
end
if numel(times) == 1 && times < 1
  permN = randperm(N);
  subN = ceil(times * N);
  times = permN(1:subN);
  times = sort(times);
end
if all(size(stations) == [2 2])
  minlon = stations(1,1);
  maxlon = stations(1,2);
  minlat = stations(2,1);
  maxlat = stations(2,2);
  stations = find(data.coordinates(1,:) >= minlon & ...
                  data.coordinates(1,:) <= maxlon & ...
                  data.coordinates(2,:) >= minlat & ...
                  data.coordinates(2,:) <= maxlat);
end

data.observations = data.observations(stations, times);
data.coordinates = data.coordinates(:,stations);
data.time = data.time(times);