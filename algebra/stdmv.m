function s = stdmv(X, dim);
% handles missing values (NaN).

error(nargchk(1,2,nargin));

if nargin < 2
  dim = min(find(size(X)~=1));
  if isempty(dim), dim = 1; end
end

sz = ones(length(size(X)), 1);
sz(dim) = size(X,dim);

% NaNs to zero
Imv = isnan(X);
X0 = X - repmat(meanmv(X,dim), sz);
X0(Imv) = 0;

% Sum and divide by the number of non-missing values
s = sum(X0.^2, dim) ./ sum(~Imv, dim);
%mu(isinf(xmean)) = NaN;
