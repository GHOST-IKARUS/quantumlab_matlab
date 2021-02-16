% ex  Expectation value of an operator
%    ex(op,rho) is the expectation value of the operator op
%    for the density matrix rho.  If rho is not normalized,
%    it is automatically normalized to have trace 1.
%    If instead of rho a state vector is given then it is 
%    automatically converted to a normalized density matrix.

function e=ex(op,rho)
   if min(size(rho))==1,
       % Input is a state vector
       e=bra(rho)*op*ket(rho);
   else
       % Input is a density matrix
       e=trace(rho*op)/trace(rho);
   end %if