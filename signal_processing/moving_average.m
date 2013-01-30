% Y = MOVING_AVERAGE(X, L)
%
% X is the original signal (column vector).
% L is the number of preceding/succeeding elements to take into account.
%
% Thus, the average is evaluated over 2xL+1 elements of X.
%
% NaNs are handled properly. If there are more than 2xL+1 consecutive NaNs,
% no mean for the corresponding element can be evaluated, and NaN is the
% result. To completely avoid NaNs in the result, use L such that 2xL+1 is
% larger than the maximum number of consecutive NaNs.
%
% Jaakko Luttinen 2009
function y = moving_average(x, l)

% Width of the filter window
d = 2*l + 1;
w = ones(d,1);

% Ignore nans
b = ~isnan(x);
x(~b) = 0;

% Add stuff to avoid bad tail behaviour
z = zeros(l,1);
x = [x(:); z];
b = [b(:); z];

% Filter, and take only the relevant stuff
y = filter(w,1,x) ./ filter(w,1,b);
y = y((l+1):end);