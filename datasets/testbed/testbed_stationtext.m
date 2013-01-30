function testbed_stationtext(Q, varargin)

options = struct(...
    'stations', 1:size(Q.W,1), ...
    'string', {{}}, ...
    'color', 'k', ...
    'verticalalignment', 'middle');

% Check arguments
[options, errmsg, remopts] = argparse(options, varargin{:});
error(errmsg);

% By default, put the stations number
if isempty(options.string)
  options.string = mat2str(options.stations);
end

maptext(Q.options.user_data.coordinates(options.stations,1:2)', options.string, ...
        'Color', options.color, 'VerticalAlignment', options.verticalalignment, ...
        remopts{:});

