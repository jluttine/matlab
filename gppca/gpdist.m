%
% D = gpdist(x1, x2, funcDist)
%
% Returns matrix of pair-wise distances between x1 and x2 as a matrix.
% x1 is K x M matrix
% x2 is K x N matrix
% where K is the dimensionality, M and N are the number of samples.
% The resulting matrix D is M x N matrix.
% If no function is given for distance calculation, 2-norm is used.
function D = gpdist(x1, x2, funcDist)

if nargin < 2
  x2 = x1;
end
if isvector(x1)
  x1 = x1(:)';
  x2 = x2(:)';
end
if nargin < 3
  funcDist = @(z1,z2) norm(z1-z2);
%  funcDist = @(z1,z2) abs(z1-z2);
end

n1 = size(x1,2);
n2 = size(x2,2);
D = zeros(n1,n2);
for i=1:n1
  for j=1:n2
    D(i,j) = funcDist(x1(:,i),x2(:,j));
  end
end
