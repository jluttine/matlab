
function clim_positive(hax)

if nargin < 1
  hax = gca;
end

cl = get(hax, 'clim');
cl = max(cl);
set(hax, 'clim', [0 cl]);
