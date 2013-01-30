function ylim_scale(hax, scale)

for n=1:numel(hax)
  ylim = get(hax(n), 'ylim');
  set(hax(n), 'ylim', scale*ylim);
end