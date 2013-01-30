% MAPTEXT   Text annotations on maps.
%
%    MAPTEXT(C, 'string') adds the text in the quotes to location
%    specified by the coordinates C. C(1,:) specifies the longitudes and
%    C(2,:) the latitudes.
%
%    MAPTEXT(C, S) adds the texts in the cell array S to locations specified
%    by the coordinates C. C(1,:) specifies the longitudes and C(2,:) the
%    latitudes.
%
%    MAPTEXT(..., property/value pairs) can be used to change, e.g.,
%    fontsize, weight, and color using the standard TEXT properties.
%
%    See also TEXT and M_TEXT.

% Copyright (c) 2010 Jaakko Luttinen

function maptext(C, str, varargin)

if iscell(str)
  for i=1:numel(str)
    m_text(C(1,i), C(2,i), str{i}, varargin{:});
  end
else
  m_text(C(1,:), C(2,:), str, varargin{:});
end
