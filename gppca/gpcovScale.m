
% K = gpcovConstScale(s, covfunc, ...)
% [K, dK_logtheta] = gpcovConstScale(s, covfunc, ...)
% [K, dK_logtheta, dK_x2] = gpcovConstScale(s, covfunc, ...)
%
% For instance,
% covfunc = {@gpcovConstScale, 0.5, {@gpcovSEiso}}
% and then
% K = feval(covfunc{:}, logtheta, x1, x2)

function varargout = gpcovScale(covfunc,logtheta,varargin)
% $$$ function [K, dK_logtheta, dK_x2] = gpcovConstScale(s,covfunc,varargin)

if ~iscell(covfunc), covfunc = {covfunc}; end

varargout = cell(nargout,1);
s = exp(logtheta(1));
logtheta(1) = []; % Don't use 2:end, because it might be empty?
[varargout{:}] = feval(covfunc{:}, logtheta, varargin{:});

for i=1:nargout
  varargout{i} = s * varargout{i};
end

% Add the partial derivative w.r.t. the scale
if nargout >= 2
  varargout{2} = cat(3, full(varargout{1}), varargout{2});
end

% $$$ figure
% $$$ in_gpcovScale = true
% $$$ imagesc(varargout{1})

% $$$ % Use constant scale
% $$$ if nargout == 1
% $$$   K = feval(covfunc{:}, varargin{:});
% $$$ elseif nargout == 2
% $$$   [K, dK_logtheta] = feval(covfunc{:}, varargin{:});
% $$$   dK_logtheta = s*dK_logtheta;
% $$$ elseif nargout == 3
% $$$   [K, dK_logtheta, dK_x2] = feval(covfunc{:}, varargin{:});
% $$$   dK_logtheta = s*dK_logtheta;
% $$$   dK_x2 = s*dK_x2;
% $$$ end
% $$$ 
% $$$ K = s * K;
