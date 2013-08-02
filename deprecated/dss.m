function [Af,Sf,S,Sv,Vf] = dss(A, S, Sv, cutoff, nComp)
% [Af,Sf,S,Vf] = dss(A, S, Sv, filter, nComp)

% $$$ error(nargchk(4,5,nargin));
% $$$ if nargin == 4
% $$$   filter = Sv;
% $$$   nComp = filter;
% $$$   Sv = [];
% $$$ end

% Denoise
%filter = filter(:)';
n = size(S,2);
[b] = fir1(floor(n/4), cutoff);
%figure
%plot(b)
Sf = filtfilt(b,1,S')';

% Identify directions with PCA
Vf = princomp(Sf');
Vf=Vf(:,1:nComp);

% Rotate
Af = A * Vf;
Sf = Vf' * Sf;
S = Vf' * S;
for j = 1:length(Sv)
    Sv{j} = Vf'*Sv{j}*Vf;
end

