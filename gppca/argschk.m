function [ opts, errmsg, wrnmsg ] = argschk( defopts, varargin )

errmsg = '';
wrnmsg = '';
opts = defopts;

if mod( length(varargin), 2 )
    if isstruct(varargin{1})
        opts_in = varargin{1};
    else
        errmsg = [ 'The optional parameters must be given in',...
                   'parameter/value pairs or using a structure argument.' ];
        return
    end
else
    opts_in = struct( varargin{:} );
end

flds = fieldnames( opts_in );
for i = 1:length(flds)
    if ~strcmp( flds{i}, lower(flds{i}) )
        opts_in.(lower(flds{i})) = opts_in.(flds{i});
        opts_in = rmfield( opts_in, flds{i} );
        flds{i} = lower(flds{i});
    end
end
for i = 1:length(flds)
    if ~isfield( opts, flds{i} )
        wrnmsg = [ 'Unknown parameter ''' flds{i} '''' ];
    end
    opts.(flds{i}) = opts_in.(flds{i});
end
