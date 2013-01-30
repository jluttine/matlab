% DATA = TESTBED_REMOVE_STATIONS(DATA, STATIONS)
%
% Removes stations from the data struct.
%
% DATA is a struct of Testbed data.
% STATIONS is the indeces of the stations to be removed.

% Last modified 2010-06-02
% Copyright (c) Jaakko Luttinen

function data_struct = testbed_remove_stations(data_struct, stations)

data_struct.observations(stations,:) = [];
data_struct.stations(stations) = [];
data_struct.coordinates(stations,:) = [];
