% dL = lchol_grad(dK, L)
%
% TODO: Consider writing this as a MEX-file!!!

function dL = lchol_grad_old(dK, L)

return

N = size(L,1);
if issparse(L)
  
  %[I,J,V] = find(L);
  %dL = tril(dK);
  
  I = find(L);
  dV = full(dK(I)); % ???
  
  
  ind = 1;
  for k=1:N
    % Pivot
    dV(ind) =  0.5 * (dV(ind)/L(k,k));
    ind_pivot = ind;
    ind = ind + 1;
% $$$     dL(k,k) =  0.5 * (dL(k,k)/L(k,k));
    
    % Lead column
    while ind <= length(I) && I(ind) <= k*N
      j = I(ind) - (k-1)*N; % row
      dV(ind) = (dV(ind) - L(j,k)*dV(ind_pivot)) / L(k,k);
      ind = ind + 1;
    end
% $$$     for j=(k+1):N
% $$$       dL(j,k) = (dL(j,k) - L(j,k)*dL(k,k)) / L(k,k);
% $$$     end
    
    % Rows
    jnd = ind;
    while jnd <= length(I)
      %dV(jnd) = 
    end
% $$$     for j=(k+1):N
% $$$       for i=j:N
% $$$         dL(i,j) = dL(i,j) - dL(i,k)*L(j,k) - L(i,k)*dL(j,k);
% $$$       end
% $$$     end

% $$$     for j=(k+1):N
% $$$       for i=j:N
% $$$         dL(i,j) = dL(i,j) - dL(i,k)*L(j,k) - L(i,k)*dL(j,k);
% $$$       end
% $$$     end
    
  end
  
  i = mod(I-1,N)+1;
  j = ceil(I/N);
  dL = sparse(i,j,dV,N,N);
  
% $$$   U = L';
% $$$   dU = spalloc(N,N,nnz(U));
% $$$   
% $$$   dU(1,1) = 0.5 * dK(1,1) / U(1,1);
% $$$   
% $$$   for j=2:N
% $$$     %dL = dU';
% $$$     i = j-1;
% $$$     z = dK(:,j) - dU'*U(:,j);
% $$$     dU(1:i,j) = linsolve_tril(L(1:i,1:i), z(1:i));
% $$$     dU(j,j) = (0.5*dK(j,j) - U(1:i,j)'*dU(1:i,j)) / L(j,j);
% $$$   end
% $$$   dL = dU';
  
else
  dU = zeros(N,N);
  
  error('not yet implemented for full matrices')
end


% $$$ for j=2:N
% $$$   i = j-1;
% $$$   z = dK(j,1:i)' - dL(1:i,1:i)*L(j,1:i)';
% $$$   dL(j,1:i) = linsolve_tril(L(1:i,1:i), z);
% $$$   dL(j,j) = (0.5*dK(j,j) - L(j,1:i)*dL(j,1:i)') / L(j,j);
% $$$ end