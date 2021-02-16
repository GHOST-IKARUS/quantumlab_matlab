%comm Commutator of two matrices.
%   comm(A,B)=AB-BA gives the commutator of A and B.

function c=comm(A,B)
c=A*B-B*A;