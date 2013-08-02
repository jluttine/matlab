
function test_argparse
warning('This function is deprecated')


options = struct('Eka', 1, 'Toka', 2);

% Simple test
[opts, errmsg] = argparse(options, 'Toka', 3)

% Test case sensitivity
[opts, errmsg] = argparse(options, 'toka', 4)

% Test unknown option with remainder
[opts, errmsg, rem] = argparse(options, 'Uusi', 10)

% Test unknown option without remainder
[opts, errmsg] = argparse(options, 'Uusi', 10)

myplot(1:10, (1:10).^2, 'Color', 'r', 'Title', 'Parabola');

myfun('Age', 99);

function myplot(x, y, varargin)
opts = struct('title', 'My magnificent plot');
[opts, errmsg, remopts] = argparse(opts, varargin{:});
error(errmsg);
plot(x, y, remopts{:});  % pass the remaining parameters to 'plot'
title(opts.title);

function myfun(varargin)
opts = struct('name', 'Foo', 'age', 42);    % default values
[opts, errmsg] = argparse(opts, varargin{:}); % parse custom values
error(errmsg);                              % print errors
fprintf('%s is %d years old.\n', opts.name, opts.age);
