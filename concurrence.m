% concurrence    Compute the concurrence of a two-qubit density matrix.

function c=concurrence(rho);
  y=[0 -i; +i 0];          
  R=rho*kron(y,y)*conj(rho)*kron(y,y);
  % Real is needed since MATLAB always adds a small imaginary part ...
  e=real(sqrt(eig(R)));
  e=-sort(-e);
  c=max(e(1)-e(2)-e(3)-e(4),0);