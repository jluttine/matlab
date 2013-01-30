function ylim_centered(hax)

for n=1:numel(hax)
  ylim = get(hax(n), 'ylim');
  ylim = max(abs(ylim));
  set(hax(n), 'ylim', [-ylim ylim]);
end