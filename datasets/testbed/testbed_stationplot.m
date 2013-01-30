function testbed_stationplot(Q, varargin)

options = struct(...
    'stations', 1:size(Q.W,1), ...
    'linespec', 'x');

% Check arguments
[options, errmsg, remopts] = argparse(options, varargin{:});
error(errmsg);

mapplot(Q.options.user_data.coordinates(options.stations,1:2)', options.linespec, ...
        remopts{:});

