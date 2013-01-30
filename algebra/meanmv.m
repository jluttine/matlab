function mu = meanmv(X, dim);
% handles missing values (NaN).

error(nargchk(1,2,nargin));

if nargin < 2
  dim = min(find(size(X)~=1));
  if isempty(dim), dim = 1; end
end

% NaNs to zero
Imv = isnan(X);
X(Imv) = 0;

% Sum and divide by the number of non-missing values
div = sum(~Imv, dim);
div(div==0) = 1;
mu = sum(X, dim) ./ div;
%mu(isinf(xmean)) = NaN;
