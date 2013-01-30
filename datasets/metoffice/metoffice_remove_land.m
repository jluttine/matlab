function X = metoffice_remove_land(X)

X(~metoffice_get_mask(),:) = [];
%S = load('/share/climate/data/UK_Met_Office/RecTest/land_hadgem');
%land = (S.land(:) >= 0.5);
%X(land,:) = [];
