% KRONPROD - Evaluates a matrix-vector product when the matrix is a
%            Kronecker product of two matrices.
%
% Y = KRONPROD(A,B,X)
%
% Evaluates KRON(A,B)*VEC(X), where VEC(X)=X(:), as Y=B*X*A'.

%
% Y = KRONPROD(...,FORM)
%
% FORM must be either 'vector' or 'matrix'. It specifies in which form the
% resulting Y is given. The default is 'matrix'.

% Last modified 2010-11-11
% Copyright (c) Jaakko Luttinen (jaakko.luttinen@tkk.fi)

function Y = kronprod(A, B, X)
% $$$ function Y = kronprod(A, B, X, form)

%X(:,:) = B*X;
%X(:,:) = X*A';

%X(:,:) = B*X*A';
Y = B*X*A';

return


% $$$ if isnumeric(B)
% $$$   Y = B*X;
% $$$ else
% $$$   Y = B(X);
% $$$ end
% $$$ if isnumeric(A)
% $$$   Y = Y*A';
% $$$ else
% $$$   Y = A(Y')';
% $$$ end
% $$$ 
% $$$ % Y = B*X*A'; % original
% $$$ 
% $$$ if nargin >= 4 && ~isempty(form)
% $$$   if strcmpi(form,'vector')
% $$$     Y = Y(:);
% $$$   elseif ~strcmpi(form,'matrix')
% $$$     error(['Requested the result in unknown form. Must be either ''vector'' ' ...
% $$$            'or ''matrix''']);
% $$$   end
% $$$ end
% $$$ 
% $$$ 
% $$$   
