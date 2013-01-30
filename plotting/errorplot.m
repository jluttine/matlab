
% ERRORPLOT(Y, E)
%
% ERRORPLOT(X, Y, E)
%
% ERRORPLOT(X, Y, L, U)
%
% X is a vector of inputs (e.g., time).
% Y is a vector of function values.
% E is a vector or a scalar defining the error as [Y-E, Y+E].
% L and U are vectors or scalars defining the error range as [Y-L, Y+U].
%
% Note that the function can plot only one signal. 
%
% Optional arguments for controlling the appearance of the error
%
%  'ErrorEdgeColor' : default is [0.8 0.8 0.8]
%  'ErrorFillColor' : default is [0.8 0.8 0.8]
%  'ErrorLineWidth' : default is 1
%  'ErrorLineStyle' : default is '-'
%  'ErrorFillAlpha' : default is 1
%  'ErrorEdgeAlpha' : default is 1
%
% In addition, the optional arguments of PLOT can be used to control the
% appearance of the function values.

% Copyright (c) 2010 Jaakko Luttinen

function hax = errorplot(y, e, varargin)

% Check whether the input was (Y,E) or (X,Y,E) or (X,Y,L,U)
if nargin >= 3 && isnumeric(varargin{1})
  x = y;
  y = e;
  if nargin >= 4 && isnumeric(varargin{2})
    % (X,Y,L,U)
    l = varargin{1};
    u = varargin{2};
    params = varargin(3:end);
  else
    % (X,Y,E)
    l = varargin{1};
    u = varargin{1};
    params = varargin(2:end);
  end
else
  % (Y,E)
  x = 1:length(y);
  l = e;
  u = e;
  params = varargin;
end

% Check that X is valid
if isvector(x)
  x = x(:);
else
  error('Input X must be a vector.');
end

% Check that Y is valid
if isvector(y)
  if length(y) == length(x)
    y = y(:);
  else
    error('Inputs X and Y must have the same length.');
  end
else
  error('Input Y must be a vector.');
end

% Check that E is valid
if isvector(l) && isvector(u)
  l = l(:);
  u = u(:);
else
  error('Errors must be given as vectors.')
end

% Parse optional arguments
options = struct(...
    'Color',          [0 0 0],          ...
    'ErrorEdgeColor', [0.8 0.8 0.8], ...
    'ErrorFillColor', [0.8 0.8 0.8], ...
    'ErrorLineWidth', 1,                ...
    'ErrorLineStyle', '-',              ...
    'ErrorFillAlpha', 1,                ...
    'ErrorEdgeAlpha', 1);
[options, errmsg, remopts] = argparse(options, params{:});
error(errmsg);

% Plot error area
X = [x(1:end); x(end:-1:1)];
Y = [y(1:end)+u; y(end:-1:1)-l(end:-1:1)];
fill(X, Y, options.ErrorFillColor, ...
     'EdgeColor', options.ErrorEdgeColor, ...
     'LineWidth', options.ErrorLineWidth, ...
     'LineStyle', options.ErrorLineStyle, ...
     'FaceAlpha', options.ErrorFillAlpha, ...
     'EdgeAlpha', options.ErrorEdgeAlpha);

% Plot the function
hold on
hax = plot(x, y, 'Color', options.Color, remopts{:});

if nargout < 1
  clear hax;
end
