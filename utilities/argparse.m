% ARGPARSE   Parse parameter and value pairs.
%
%    [OPTIONS, ERROR_MSG] = ARGPARSE(DEFAULT_OPTIONS, ...)
%
%    DEFAULT_OPTIONS is a struct of parameters and their default values. The
%    custom parameter values are given as a struct or as parameter/value
%    pairs.
%
%    ERROR_MSG contains an parsing error message as a string or is empty if
%    no errors occured.
%
%    For example,
%
%      function myfun(varargin)
%      opts = struct('name', 'Foo', 'age', 42);      % default values
%      [opts, errmsg] = argparse(opts, varargin{:}); % parse custom values
%      error(errmsg);                                % print errors
%      fprintf('%s is %d years old.\n', opts.name, opts.age);
%
%    Now,
%
%      myfun('Age', 99)
%
%    results
%
%      Foo is 99 years old.
%
%    By default, the function returns an error message if the list of custom
%    parameter values contains parameter names not included in the list of
%    default options. This can be prevented by using a third output argument
%    as:
%
%    [OPTIONS, ERROR_MSG, REMAINING_OPTIONS] = ARGPARSE(DEFAULT_OPTIONS, ...)
%
%    REMAINING_OPTIONS is a cell array of parameters and their values for
%    such parameters that were not contained in DEFAULT_OPTIONS. This is
%    useful for forwarding the remaining unknown arguments to some other
%    function.
%
%    For example,
%
%      function myplot(x, y, varargin)
%      opts = struct('title', 'My magnificent plot');
%      [opts, errmsg, remopts] = argparse(opts, varargin{:});
%      error(errmsg);
%      plot(x, y, remopts{:});  % pass the remaining parameters to 'plot'
%      title(opts.title);
%
%    Now,
%
%      myplot(1:10, (1:10).^2, 'Color', 'r', 'Title', 'Parabola');
%
%    would result a plot of a parabola with red color and a title
%    'Parabola'.
%
%    Note: The parsing of arguments is case insensitive. Thus,
%    DEFAULT_OPTIONS should not contain such field names that differ only in
%    case.

% Last modified 2011-10-20
% Copyright (c) Jaakko Luttinen (jaakko.luttinen@aalto.fi)
%
% Based on code by Alexander Ilin and Tapani Raiko.

function [opts, error_msg, rem_opts] = argparse( defopts, varargin )

error( nargoutchk(2,3,nargout) );

error_msg = '';
rem_opts = {};

opts = defopts;

if nargin == 2 && isempty(varargin{1})
  return
end

% Convert argument pairs to struct
if mod( length(varargin), 2 )
  if isstruct(varargin{1})
    opts_in = varargin{1};
  else
    error_msg = ['The optional parameters must be given in parameter/value' ...
                 ' pairs or using a structure argument.'];
    return
  end
else
  opts_in = struct( varargin{:} );
end

% Set matching options
fields_in = fieldnames(opts_in);
for i = 1:length(fields_in)
  fieldname = matchfield(opts, fields_in{i});
  value = getfield(opts_in, fields_in{i});
  if ~isempty(fieldname)
    if ischar(fieldname)
      % Found a matching field
      opts.(fieldname) = value;
      opts = setfield(opts, fieldname, value);
    else
      % Several fields matched
      error_msg = ['Could not uniquely resolve argument ''' fields_in{i} ''''];
    end
  else
    % No matching fields
    if nargout < 3
      error_msg = [ 'Unknown parameter ''' fields_in{i} '''' ];
    else
      rem_opts = {rem_opts{:}, fields_in{i}, value};
    end
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function field = matchfield(s, name)

exact_match = [];
partial_match = [];

L = length(name);

names = fieldnames(s);
for n = 1:length(names)
  if strcmpi(names{n}, name)
    if isempty(exact_match)
      % Unique perfect match
      exact_match = names{n};
    else
      % Error: Two matching fields
      field = -1;
      return
    end
  elseif strncmpi(names{n}, name, L)
% $$$     if isempty(partial_match)
% $$$       % Unique partial match
% $$$       partial_match = names{n};
% $$$     else
% $$$       % Error: Two partially matching fields
% $$$       field = -1;
% $$$       return
% $$$     end
  end
end

field = [];
if ~isempty(exact_match)
  field = exact_match;
elseif ~isempty(partial_match)
  field = partial_match;
end
