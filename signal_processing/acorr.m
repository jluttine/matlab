function c = acorr(y, maxlag)

if isvector(y)
  if nargin < 2
    maxlag = length(y) - 1;
  end
  y = y - mean(y);
  c = zeros(1+maxlag,1);
  c0 = acov(y,0)
  c(1) = 1;
  for l=1:maxlag
    c(l+1) = acov(y, l) / c0;
  end
else
  if nargin < 2
    maxlag = size(y,1) - 1;
  end
  y = bsxfun(@minus, y, mean(y,1));
  c = zeros(maxlag+1,size(y,2));
  for n=1:size(y,2)
    c0 = acov(y(:,n),0);
    c(1,n) = 1;
    for l=1:maxlag
      c(l+1,n) = acov(y(:,n),l) / c0;
    end
  end
end

function c = acov(x0, L)
N = length(x0);
c = 1/(N-L) * ( dot(x0(1:(N-L)), x0((1+L):end)) );
