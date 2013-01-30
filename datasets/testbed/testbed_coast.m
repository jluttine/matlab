% TESTBED_COAST - Draws the coast line of the Testbed research area.
%
% TESTBED_COAST()
% TESTBED_COAST(HAX)
%
% HAX : a vector of axes handles [default: GCA]

% Last modified 2010-06-11
% Copyright (c) Jaakko Luttinen (jaakko.luttinen@tkk.fi)

function testbed_coast(hax)

if nargin < 1
  hax = gca;
end

map_coast(hax, 'usercoast','testbed_data_coast');
