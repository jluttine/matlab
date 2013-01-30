% H = XLABELS(TEXT,...)
% H = XLABELS(HAX,TEXT,...)

% Last modified 2010-06-11
% Copyright (c) Jaakko Luttinen (jaakko.luttinen@tkk.fi)

function h = xlabels(text,varargin)

options = struct('Separation', []);

if isnumeric(text)
  hax = text;
  text = varargin{1};
  [options, errmsg, remopts] = argparse(options, varargin{2:end});
  error(errmsg);
else
  hax = gca;
  [options, errmsg, remopts] = argparse(options, varargin{:});
  error(errmsg);
end

h = zeros(size(hax));
for n=1:numel(hax)

  % Get axes position
  pos_ax = get(hax(n), 'Position');
  
  % Add a colorbar
  if iscell(text)
    h(n) = xlabel(hax(n), text{n}, remopts{:});
  else
    h(n) = xlabel(hax(n), text, remopts{:});
  end
  
  % Use the same units
  units = get(hax(n), 'Units');
  set(h(n), 'Units', units);
  
  pos = get(h(n), 'Position');
  pos(2) = -options.Separation;
  set(h(n), 'Position', pos);

end

if nargout < 1
  clear hcb;
end