function set_plot_style(hax)
if nargin < 1
  hax = gca;
end
set(hax, 'Box', 'off' ); 
set(hax, 'TickDir', 'out'); 