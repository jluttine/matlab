% K = gpcovProduct(covfunc1, covfunc2, ...)
% [K, dK_logtheta] = gpcovConstScale(covfunc1, covfunc2, ...)
% [K, dK_logtheta, dK_x2] = gpcovConstScale(covfunc1, covfunc2, ...)
%
% For instance,
% covfunc = {@gpcovProduct, {@gpcovSEiso}, {@gpcovPeriodic}}
% and then
% K = feval(covfunc{:}, logtheta, x1, x2)

function varargout = gpcovProduct(covfunc1,covfunc2,logtheta,varargin)
% $$$ function [K, dK_logtheta, dK_x2] = gpcovConstScale(s,covfunc,varargin)

if ~iscell(covfunc1), covfunc1 = {covfunc1}; end
if ~iscell(covfunc2), covfunc2 = {covfunc2}; end

params1 = feval(covfunc1{:});
params2 = feval(covfunc2{:});

logtheta1 = logtheta(1:params1);
logtheta2 = logtheta((params1+1):end);

varargout1 = cell(nargout,1);
[varargout1{:}] = feval(covfunc1{:}, logtheta1, varargin{:});
varargout2 = cell(nargout,1);
[varargout2{:}] = feval(covfunc2{:}, logtheta2, varargin{:});

varargout = cell(nargout,1);
if nargout >= 1
  varargout{1} = varargout1{1}.*varargout2{1};
end
if nargout >= 2
  p = params1 + params2;
  [m,n] = size(varargout{1});
  varargout{2} = zeros([m,n,p]);
  varargout{2}(:,:,1:params1) = bsxfun(@times,varargout1{2},varargout2{1});
  varargout{2}(:,:,(params1+1):end) = bsxfun(@times,varargout1{1},varargout2{2});
end
if nargout >= 3
  [m,n] = size(varargout{1});
  varargout{3} = bsxfun(@times,reshape(varargout1{1},[1,m,n]),varargout2{3}) + ...
      bsxfun(@times,varargout1{3},reshape(varargout2{1},[1,m,n]));
end

% $$$ % Add the partial derivative w.r.t. the scale
% $$$ if nargout >= 2
% $$$   varargout{2} = cat(3, full(varargout{1}), varargout{2});
% $$$ end

