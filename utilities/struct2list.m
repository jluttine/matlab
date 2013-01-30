% L = STRUCT2LIST(S)
%
% Converts a struct into list of fieldname/value pairs.
%
% For instance,
%
%   s = struct('name', 'Foo Bar', 'age', 42);
%
% Then
% 
%   l = struct2list(s);
%
% is equivalent to
%
%   l = {'name', 'Foo Bar', 'age', 42};
%
% Thus,
%
%   struct(struct2list(s))
%
% is identity operator.

function l = struct2list(s)

names = fieldnames(s);
values = struct2cell(s);

N = length(names);
l = cell(2*N,1);

for n=1:N
  i = 2*n-1;
  l(i:(i+1)) = {names{n}, values{n}};
end