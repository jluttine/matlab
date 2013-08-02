function X = solve_tril(U, Y)
% function X = solve_tril(L, Y)
%
% Solves X = L \ Y, or L*X = Y,  when L is lower triangular.

if issparse(U)
  X = U \ Y;
else
  opts.UT = true;
  X = linsolve(U, full(Y), opts);
end
