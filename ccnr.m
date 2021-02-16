% ccnr   Computable cross norm - realignment criterion
%    ccnr(rho) gives the trace norm of the realigned rho.
%    If it is larger than 1 then the system is entangled.
%    rho must be bipartite, with two parties of equal dimensions.

function c=ccnr(rho)

rho=ketbra2(rho);
[sy,sx]=size(rho);
c=trnorm(realign(rho));



