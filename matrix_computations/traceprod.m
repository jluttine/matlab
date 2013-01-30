% X = TRACEPROD(A, B)
%
% Evaluates trace(A*B') efficiently.

% Last modified 2010-06-08
% Copyright (c) Jaakko Luttinen (jaakko.luttinen@tkk.fi)

function x = traceprod(A,B)


if any(size(A)~=size(B))
  error('The input matrices must have the same number of rows and columns.');
end

if nargin < 2
  warning('Do not use like this.');
  x = A(:)'*A(:);
else
  if issparse(A) || issparse(B)
    x = full(sum(sum(A.*B)));
  else
    x = A(:)'*B(:);
  end
end