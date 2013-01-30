
function L = ldlchol2lchol(LD)
[L,D] = ldlsplit(LD);
L = L * sqrt(D);
