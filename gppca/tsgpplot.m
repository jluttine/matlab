
function hax = tsgpplot(x, Y, E, varargin)
% function tsgpplot(x, Mu, E)
%
% Plots d signals.
% x = inputs (n x 1)
% Y = expected outputs (n x d)
% E = errors of the outputs (n x d)

opts = struct('pseudoinputs', []);
[ opts, errmsg, wrnmsg ] = argschk( opts, varargin{:} );
if ~isempty(errmsg), error( errmsg ), end
if ~isempty(wrnmsg), warning( wrnmsg ), end

figure

n = size(Y,2);
for i=1:n
  params = {};
  if ~isempty(opts.pseudoinputs)
    params = [params(:); {'pseudoinputs'}; {opts.pseudoinputs{i}}];
  end
  hax(i) = subplot(n,1,i);
  gpplot(x(:),Y(:,i),E(:,i), params{:});
  
  set( gca, 'xlim', [ min(x) max(x) ] )
  
  pos_ax = get( gca, 'Position' );
  hgt_ax = (0.95 - 0.1) / ( n + (n-1)*0.1 );
  hgt_sp = hgt_ax * 0.1;
  pos_ax(1) = 0.1;
  pos_ax(2) = 0.95 - (i-1)*( hgt_ax + hgt_sp ) - hgt_ax;
  pos_ax(3) = 0.84;
  pos_ax(4) = hgt_ax;
  set( gca, 'Position', pos_ax )
  if i ~= n
    disp('removing x tick labels')
    set( gca, 'XTickLabel', [] )
  end
  axis tight

end

if nargout < 1
  clear hax;
end
  