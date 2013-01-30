% LOGDET_DIAG - Computes the log-determinant of a diagonal matrix.
%
% Y = LOGDET_DIAG(D)
%
% D must be a diagonal matrix.

% Last modified 2011-01-24
% Copyright (c) Jaakko Luttinen

function y = logdet_diag(D)

y = sum(log(diag(D)));