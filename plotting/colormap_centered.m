function colormap_centered(hax)

if nargin < 1
  hax = gca();
end

for n=1:length(hax)
  cl = get(hax(n), 'clim');
  l = max(abs(cl));
  cl = [-l l];
  set(hax(n), 'clim', cl);
end