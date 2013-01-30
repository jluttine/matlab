function c = cellsum(c1, c2)

c = cell(size(c1));
for n=1:numel(c1)
  c{n} = c1{n} + c2{n};
end
