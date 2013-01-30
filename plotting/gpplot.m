
function gpplot(x, y, e, varargin)
% function gpplot(x, y, e, flipaxis)
%
% x = inputs (e.g., time)
% y = expectated values of the outputs
% e = errors of the outputs

warning('Deprecated, use errorplot')

opts = struct('pseudoinputs', []);
[ opts, errmsg, wrnmsg ] = argschk( opts, varargin{:} );
if ~isempty(errmsg), error( errmsg ), end
if ~isempty(wrnmsg), warning( wrnmsg ), end



% $$$ if nargin < 4
% $$$   flipaxis = false;
% $$$ end
x = x(:);
y = y(:);
%e = sqrt(diag(Cov));
X = [x(1:end); x(end:-1:1)];
Y = [y(1:end)+e; y(end:-1:1)-e(end:-1:1)];
C = [.75,.75,.75];
% $$$ if flipaxis
% $$$   fill(Y,X,C,'EdgeColor',C);
% $$$   hold on
% $$$   plot(y,x,'k');
% $$$ else
fill(X,Y,C,'EdgeColor',C);
hold on
plot(x,y,'k');
if ~isempty(opts.pseudoinputs)
  yl = ylim;
  py = yl(1) + 0.8*(yl(2)-yl(1));
  plot(opts.pseudoinputs, py, 'k+')
end
% $$$ end