function X = metoffice_remove_bins(X,maskfile)

X(~metoffice_get_mask(maskfile),:) = [];
