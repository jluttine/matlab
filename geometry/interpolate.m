% YH = INTERPOLATE(D, Y, DH)
%
% Interpolates the values of Y based on given distances.
%
% D  : N x N matrix of distances between the given targets Y.
% Y  : N x 1 vector of function values.
% DH : M x N matrix of distances between the interpolated and given targets.
%
% Returns:
%
% YH : M x 1 vector of interpolated targets.
%
% Optional arguments:
%
% 'method' : 'nn' means nearest neighbour interpolation (default).
%            'rbf' means linear radial basis function interpolation.

% Copyright (c) 2010 Jaakko Luttinen

function yh = interpolate(D, y, Dh, varargin)

options = struct( ...
    'method', 'nn');

% Get optional arguments
[options, errmsg] = argparse(options, varargin{:});
error(errmsg);

if strcmpi(options.method, 'nn')
  
  %% Nearest neighbour interpolation
  
  [tmp,ind] = min(Dh,[],2);
  yh = y(ind);
  
elseif strcmpi(options.method, 'rbf')
  
  %% Linear radial basis functions interpolation
  
  mu = mean(y);
  y = y - mu;

  yh = Dh * (D\y) + mu;
  
else
  
  error('Unknown interpolation method requested.');
  
end
