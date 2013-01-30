function c = acorr(y, maxlag)


if isvector(y)
  if nargin < 2
    maxlag = length(y);
  end
  y = y - mean(y);
  z = xcorr(y,maxlag,'coeff');
  l0 = (length(z)+1)/2;
  c = z(l0:end);
else
  if nargin < 2
    maxlag = size(y,1);
  end
  y = bsxfun(@minus, y, mean(y,1));
  c = zeros(maxlag+1,size(y,2));
  for n=1:size(y,2)
    z = xcorr(y(:,n),maxlag,'coeff');
    l0 = (size(z,1)+1)/2;
    c(:,n) = z(l0:end);
  end
end
    