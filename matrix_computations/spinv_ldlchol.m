% SPINV_LDLCHOL - Evaluate the sparse inverse matrix given the LDL
%                 decomposition.
%
% For a sparse symmetric positive definite matrix C,
%
%   Z = SPINV(LD)
%
% where LD = LDLCHOL(C) and [Z]_ij = [inv(C)]_ij for such ij that [C]_ij
% is non-zero.
%
% See Vanhatalo and Vehtari (2008) for details. 

% Copyright (c) 2008      Jarno Vanhatalo

% This software is distributed under the GNU General Public
% License (version 2 or later); please refer to the file
% License.txt, included with the software, for details.

function Z = spinv_ldlchol(LD)

Z = spinv(LD, 1);
return

    
    n = size(LD,1);


    % The mex-file was not available, so evaluate the sparse inverse here
    
    % TODO:
    % For now, just evaluate the matrix A and then run..
    %[L,D] = ldlsplit(LD);
    %A = L*D*L';
    
    if nargin < 2
      q = 1:n;
    end

    %[LD, p, q] = ldlchol(A);
    
    [I,J,ld] = find(LD);
    temp = [I(:) J(:) ; J(:) I(:)];
    temp = sortrows(unique(temp,'rows'),2);
    Iz = temp(:,1); Jz = temp(:,2); 
    
    % Find the column starting points
    a1=zeros(n,1);
    a2 = cumsum(histc(J,1:n));
    a1(1) = 1; a1(2:end) = a2(1:end-1) + 1;
    az1=zeros(n,1);
    az2 = cumsum(histc(Jz,1:n));
    az1(1) = 1; az1(2:end) = az2(1:end-1) + 1;
    
    for j=1:n
        indaz{j} = az1(j):az2(j);
        indIz{j} = Iz(indaz{j})';
    end

    % Evaluate the sparse inverse
    z = zeros(size(Iz));
    z(end) = 1./ld(end);
    % Allocate memory
    cindit=zeros(n,1);
    for jj = n-1:-1:1
        fil = ld(a1(jj)+1:a1(jj+1)-1);
        fi = I(a1(jj)+1:a1(jj+1)-1);
        lfi = length(fi);
        Zt = zeros(lfi,lfi);
        indz = cumsum(histc(indIz{jj},[0 ; fi]));
        indz = az1(jj) + indz(1:end-1);
        
        i4=0;            
        for i1 = 1:lfi
            cind1=indaz{fi(i1)};
            Icind1=indIz{fi(i1)};
            indfi = lfi;
            i2=length(Icind1);
            go = true;
            while go
                if Icind1(i2)==jj  % Find the indeces for the jj'th rows in fi columns
                    i4=i4+1;
                    cindit(i4)=cind1(i2);
                    go = false;
                end
                if indfi >= 1 && fi(indfi) == Icind1(i2) % Find the indeces for the fi'th rows in i2'nd columns
                    Zt(indfi,i1) = z(cind1(i2));
                    indfi = indfi-1;
                end
                i2 = i2-1;
            end
        end
        % remove extras
        cindi=cindit(1:i4);

        zij = -fil'*Zt;
        z(cindi) = zij;
        z(indz) = zij;
        zij = 1./ld(a1(jj)) - fil'*z(indz);
        z(az1(jj)-1+find(indIz{jj}==jj,1)) = zij;
    end
    
    Z = sparse(Iz,Jz,z);
    r(q) = 1:n;
    Z = Z(r,r);
    
end
