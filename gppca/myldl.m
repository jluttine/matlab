function [L,D] = myldl(A)

% TODO: Check squareness and symmetricity!!


n = rows(A);

L = eye(n);
D = zeros(n);

for i=1:n
  m = i - 1;
  for j=1:m
    r = j - 1;
    L(i,j) = 1/D(j,j) * (A(i,j) - L(i,1:r)*D(1:r,1:r)*L(j,1:r)');
  end
  D(i,i) = A(i,i) - L(i,1:m)*D(1:m,1:m)*L(i,1:m)';
end
