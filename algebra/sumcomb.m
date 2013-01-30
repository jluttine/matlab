% SUMCOMB 
%
% y = sum(sum(bsxfun(@plus, x1(:), x2(:)')));

% Last modified 2010-11-09
% Copyright (c) Jaakko Luttinen (jaakko.luttinen@tkk.fi)

function y = sumcomb(x1,x2)

% there must some more efficient way..
y = sum(sum(bsxfun(@plus, x1(:), x2(:)')));
