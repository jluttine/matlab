
% DATETICKS(...)
%
% DATETICKS(AXES, ...)
%
% Simple wrapper to run DATETICK on several plots.
%
% See DATETICK for details on the optional arguments.

% Copyright (c) 2010 Jaakko Luttinen

function dateticks(varargin)

if nargin >= 1 && isnumeric(varargin{1})
  hax = varargin{1};
  params = varargin(2:end);
else
  hax = gca;
  params = varargin;
end

N = length(hax);
for n=1:N
  datetick(hax(n), params{:});
end