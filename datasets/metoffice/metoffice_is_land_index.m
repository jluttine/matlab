function b = metoffice_is_land_index(ind)

error('don''t use this function')
%data = mohsst5_loaddata();

%land = sum(~isnan(data.observations),2) == 0;
m = metoffice_get_mask();
b = ~m(ind);