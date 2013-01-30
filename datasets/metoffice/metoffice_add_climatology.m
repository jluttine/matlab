% Add climatology to the full fields
function Yrec = add_climatology( Yrec, dataset, anomalies )

if anomalies
  datapath = metoffice_get_path();
  load( [ datapath, dataset ], 'time', 'clim' )
    
  dv = datevec( time ); month = dv(:,2);
    
  Yrec = Yrec + clim(:,month);

end
