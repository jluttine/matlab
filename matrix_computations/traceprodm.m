% X = TRACEPRODM(A,B)
%
% Evaluates TRACE(A*B') efficiently.
%
% If A and/or B have a third dimension, the function returns a vector of
% traces as:
%
% X(I) = TRACE(A(:,:,I)*B(:,:,I)')
%
% A or B is kept fixed if it does not have a third dimension.

% Last modified 2010-06-08
% Copyright (c) Jaakko Luttinen (jaakko.luttinen@tkk.fi)

function x = traceprodm(A,B)

warning('Deprecated?');

[rows_A, cols_A, pages_A] = size(A);
[rows_B, cols_B, pages_B] = size(B);

if rows_A~=rows_B || cols_A~=cols_B
  error('The input matrices must have the same number of rows and columns.');
end

if pages_A~=pages_B && pages_A>1 && pages_B>1
  error('Non-singular third dimensions must agree.');
end

x = sum( bsxfun(@times, ...
                reshape(A,[rows_A*cols_A,pages_A]), ...
                reshape(B,[rows_B*cols_A,pages_B])), 1 )';
