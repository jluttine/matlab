
% TSERRORPLOT(Y, E)
%
% TSERRORPLOT(X, Y, E)
%
% TSERRORPLOT(X, Y, L, U)
%
% Plots D signals with errors.
%
% X : N x 1 matrix of inputs
% Y : N x D matrix of function values
% E : N x D or N x 1 or 1 x D or 1 x 1 matrix of errors
% L : N x D or N x 1 or 1 x D or 1 x 1 matrix of lower errors
% U : N x D or N x 1 or 1 x D or 1 x 1 matrix of upper errors
%
% See ERRORPLOT for details on the optional arguments.

% Copyright (c) 2010 Jaakko Luttinen

function hax = tserrorplot(Y, E, varargin)

% Check whether the input was (Y,E) or (X,Y,E) or (X,Y,L,U)
if nargin >= 3 && isnumeric(varargin{1})
  x = Y;
  Y = E;
  if nargin >= 4 && isnumeric(varargin{2})
    % (X,Y,L,U)
    L = varargin{1};
    U = varargin{2};
    params = varargin(3:end);
  else
    % (X,Y,E)
    L = varargin{1};
    U = varargin{1};
    params = varargin(2:end);
  end
else
  % (Y,E)
  x = (1:(size(Y,1)))';
  L = E;
  U = E;
  params = varargin;
end

% Parse optional arguments
options = struct();
[options, errmsg, remopts] = argparse(options, params{:});
error(errmsg);

% Scale L and U to proper sizes
L = bsxfun(@times, ones(size(Y)), L);
U = bsxfun(@times, ones(size(Y)), U);

figure

D = size(Y,2);
hax = zeros(D,1);

for d=1:D
  % Plot with errors
  hax(d) = subplot(D,1,d);
  errorplot(x(:),Y(:,d),L(:,d),U(:,d), remopts{:});
end

if nargout < 1
  clear hax;
end
  