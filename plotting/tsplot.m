%  TSPLOT - Plot time series
%
%  TSPLOT(Y) plots time series in the rows of Y on a separate
%  subplot. Y can be followed by parameter/value pairs to specify
%  additional properties of the lines.
%
%  See also ADDTSPLOT

%  This software is provided "as is", without warranty of any kind.
%  Alexander Ilin, Tapani Raiko

function hax = tsplot( y, varargin )

if nargin >= 2 & isnumeric(varargin{1})
    x = y;
    y = varargin{1};
    varargin = {varargin{2:end}};
else
    x = 1:size(y,2);
end

nts = size(y,1);

figure
for i = 1:nts
    hax(i) = subplot( nts, 1, i );
end

for i = 1:nts
    axes(hax(i))
    
    plot( x, y(i,:), varargin{:} )
    set( gca, 'xlim', [ min(x) max(x) ] )
    
    pos_ax = get( gca, 'Position' );
    
    hgt_ax = (0.95 - 0.1) / ( nts + (nts-1)*0.1 );
    hgt_sp = hgt_ax * 0.1;
    pos_ax(1) = 0.1;
    pos_ax(2) = 0.95 - (i-1)*( hgt_ax + hgt_sp ) - hgt_ax;
    pos_ax(3) = 0.84;
    pos_ax(4) = hgt_ax;
    set( gca, 'Position', pos_ax )
    if i ~= nts
        set( gca, 'XTickLabel', [] )
    end
end

if nargout == 0
    clear hax
end