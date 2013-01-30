%
% [D2, dD1] = sqdistEuclidean(x1, x2)
%
function [D, dD_dx2] = sqdistEuclidean(x1, x2, squared) 

if nargin < 3
  squared = true;
end

D = sq_dist(x1,x2);

if nargout >= 2
  [d,n1] = size(x1);
  n2 = cols(x2);
  dD_dx2 = 2 * bsxfun(@minus, x1, reshape(x2,[d,1,n2])) * (-1);
end

%    dD_dx2(:,:,i) = 2 * bsxfun(@minus, x1, x2(:,i)) * (-1);

warning off MATLAB:divideByZero
if squared == false
  D = sqrt(D);
  if nargout >= 2
    %[i,j] = find(D>0);
    dD_dx2 = 0.5 * bsxfun(@rdivide, dD_dx2, reshape(D,[1,n1,n2]));
    dD_dx2(isnan(dD_dx2)) = 0; % TODO, correct?
  end
  
% $$$   for i=1:n2
% $$$ % $$$ %    dD_dx2 = 2 * bsxfun(@minus, x1, reshape(x2,[d,1,n2])) * (-1);
% $$$ % $$$     dD_dx2(:,:,i) = 2 * bsxfun(@minus, x1, x2(:,i)) * (-1);
% $$$ % $$$     %error('he')
% $$$ % $$$     if squared == false
% $$$ % $$$       D = sqrt(D);
% $$$ % $$$       error('not yet')
% $$$ % $$$     end
% $$$   end
end

