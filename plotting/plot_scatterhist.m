
function plot_scatterhist(x)

N = size(x,2);

for n=1:N
  for m=1:N
    subplot(N,N,n+(m-1)*N)
    if n ~= m
      plot(x(:,n),x(:,m),'k.')
    else
      hist(x(:,n));
    end
  end
end
