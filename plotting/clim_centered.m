
function clim_centered(hax)

if nargin < 1
  hax = gca;
end

cl = get(hax, 'clim');
cl = max(abs(cl));
set(hax, 'clim', [-cl cl]);
