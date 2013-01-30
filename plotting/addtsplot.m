%  ADDTSPLOT - Add time series to current figure
%
%  ADDTSPLOT(Y) adds time series in the rows of Y to the current
%  figure which should be created by TSPLOT. Y can be followed by
%  parameter/value pairs to specify additional properties of the
%  lines.
%
%  See also TSPLOT

%  This software is provided "as is", without warranty of any kind.
%  Alexander Ilin, Tapani Raiko

function hax = addtsplot( y, varargin )

if nargin >= 2 & isnumeric(varargin{1})
    x = y;
    y = varargin{1};
    varargin = {varargin{2:end}};
else
    x = 1:size(y,2);
end

nts = size(y,1);

hax = findobj( gcf, 'type', 'axes' );
if length(hax) ~= nts
    tsplot( y, varargin{:} );
    return
end

pos = zeros( length(hax), 4 );
for i = 1:length(hax)
    pos(i,:) = get( hax(i), 'position' );
end
[tmp,I] = sort( -pos(:,2) );
hax = hax(I);

for i = 1:nts
    axes(hax(i))
    
    hold on
    plot( x, y(i,:), varargin{:} )

end

if nargout == 0
    clear hax
end