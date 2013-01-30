function set_label_fontsize(hax, size)

for i=1:length(hax)
  h_xlabel = get(hax(i),'XLabel');
  set(h_xlabel,'FontSize',size);
  h_ylabel = get(hax(i),'YLabel');
  set(h_ylabel,'FontSize',size);
end
