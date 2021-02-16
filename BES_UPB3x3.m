% BES_UPB3x3   3x3 bound entangled state constructed with UPBs
%    BES_UPB3x3 gives a two-qutrit bound entangled state with 
%    constructed with unextendible product bases.
%    See quant-ph/9908070v3.

function rho=BES_UPB3x3

f0=ket(kron([1  0  0],[1 -1  0]));
f1=ket(kron([1 -1  0],[0  0  1]));
f2=ket(kron([0  0  1],[0  1 -1]));
f3=ket(kron([0  1 -1],[1  0  0]));
f4=ket(kron([1  1  1],[1  1  1]));

rho=nm(qeye(2,3)-ketbra(f0)-ketbra(f1)-ketbra(f2)-ketbra(f3)-ketbra(f4));



