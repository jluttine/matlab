% SET_COLORBAR_POSITION
%
% SET_COLORBAR_POSITION(HCB, HAX, CB_SIZE, CB_SEP)
%
% HCB is the handle to the colorbar object
%
% HAX is the handle to the axes
%
% CB_SIZE is the size of the colorbar
%
% CB_SEP is the separation of the colorbar from the axes
%
% At the moment, this only works for horizontal colorbars placed under the
% axes.

function set_colorbar_position(hcb, hax, cb_size, cb_sep)

for d=1:numel(hax)
  pos_ax = get( hax(d), 'Position');
  cb_pos = pos_ax;
  cb_pos(2) = pos_ax(2) - cb_size - cb_sep;
  cb_pos(4) = cb_size;
  set(hcb(d), 'Position', cb_pos);
end
