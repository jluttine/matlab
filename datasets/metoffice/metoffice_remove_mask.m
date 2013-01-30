function X = metoffice_remove_mask(X,dataset)

X(~metoffice_get_mask(dataset),:) = [];

function X = metoffice_remove_bins_old(X,maskfile)

X(~metoffice_get_mask(maskfile),:) = [];
