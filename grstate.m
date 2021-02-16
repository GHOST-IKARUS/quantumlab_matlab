% grstate   Normalized ground state of a Hamiltonian

function gs=grstate(H)
   [V,D]=eig(H);
   [junk,index]=min(real(diag(D)));
   gs=V(:,index(1));
   gs=gs/sqrt(gs'*gs);