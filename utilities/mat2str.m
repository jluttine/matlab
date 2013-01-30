% MAT2STR   Convert numbers to strings.
%
%    Y = MAT2STR(X) converts the elements of the matrix X to strings. The
%    strings are contained in a cell array with the same size as X.

% Copyright (c) 2010 Jaakko Luttinen

function Y = mat2str(X)

Y = cell(size(X));

for i=1:numel(X)
  Y{i} = num2str(X(i));
end
