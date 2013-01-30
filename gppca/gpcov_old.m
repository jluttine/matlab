function [K, dK] = gpcov(varargin)

K = 0;
dK = [];
for i=1:nargin
  if nargout == 1
    Knew = varargin{i}();
  else
    [Knew,dKnew] = varargin{i}();
    n = size(dKnew,3);
    if isempty(dK)
      dK = dKnew;
    else
      dK(:,:,end+(1:n)) = dKnew;
    end
  end
  K = K + Knew;
end
