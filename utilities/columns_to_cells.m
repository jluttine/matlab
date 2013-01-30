% COLUMNS_TO_CELLS - Puts the columns of a matrix to cells
%
%   Y = COLUMNS_TO_CELLS(X)
%
% An MxN matrix is transformed to Nx1 cell array, where each cell
% contains a column vector of X: Y{K} = X(:,K).

% Last modified 2011-01-27
% Copyright (c) Jaakko Luttinen (jaakko.luttinen@tkk.fi)

function y = columns_to_cells(X)
N = size(X,2);
y = cell(N,1);
for n=1:N
  y{n} = X(:,n);
end
  