function [T,p] = statchol(Sigma,flag)
warning('This function is deprecated')

%STATCHOL  Do Cholesky-like decomposition, allowing zero eigenvalues
%   T = STATCHOL(SIGMA) computes T such that SIGMA = T'*T.  Sigma must be
%   square and symmetric.  If SIGMA is positive definite, then T is the
%   square, upper triangular Cholesky factor.  In general however, T is not
%   square or triangular.
%
%   [T,P] = STATCHOL(SIGMA) returns the number of number of negative
%   eigenvalues of SIGMA, and T is empty if P>0.  If P==0, SIGMA is positive
%   semi-definite.
%
%   If SIGMA is not square and symmetric, P is NaN and T is empty.
%
%   [T,P] = STATCHOL(SIGMA,0) returns P==0 if SIGMA is positive definite, and
%   T is the Cholesky factor.  If SIGMA is not positive definite, P is a
%   positive integer and T is empty.  [...] = STATCHOL(SIGMA,1) is equivalent
%   to [...] = STATCHOL(SIGMA).

%   Copyright 1993-2005 The MathWorks, Inc. 
%   $Revision: 1.1.4.3 $  $Date: 2005/11/18 14:28:49 $

if nargin < 2, flag = 1; end

% Test for square, symmetric
[n,m] = size(Sigma);
%%%if (n == m) & all(all(abs(Sigma - Sigma') < 10*eps(max(abs(diag(Sigma)))))) % Original IF statement
if (n == m) & all(all(abs(Sigma - Sigma') < 10*eps(full(max(abs(diag(Sigma))))))) % MODIFIED to handle sparse matrices
    [T,p] = chol(Sigma);

    % Test for positive definiteness
    if p > 0 && flag
        % Can get factors of the form Sigma==T'*T using the eigenvalue
        % decomposition of a symmetric matrix, so long as the matrix
        % is positive semi-definite.
        [U,D] = eig((Sigma+Sigma')/2);
        D = diag(D);
        
        tol = eps(max(D)) * length(D);
        t = (abs(D) > tol);
        D = D(t);
        p = sum(D<0); % number of negative eigenvalues

        if (p==0)
            T = diag(sqrt(D)) * U(:,t)';
        else
            T = [];
        end
    end

else
    T = [];
    p = NaN;
end

