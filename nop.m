% nop   Gives the anumber operator as a matrix. Usage: nop(d),
%       where d is the size of the matrix.

function m=nop(d);

m=diag(0:d-1);