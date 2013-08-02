
function [Y1, Y2] = split_data(Y, I)
warning('This function is deprecated')


Y1 = Y;
Y2 = Y;

Y1(~I) = nan;
Y2(I) = nan;


