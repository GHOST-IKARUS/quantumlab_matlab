% ising_classical_ground   Ground state energy of the classical Ising model 
%   ising_classical_ground(B) gives the ground state energy per spin for 
%   the classical Ising model in transverse field if the 
%   field strength is B and the coefficient 
%   of the nearest neighbor coupling is 1.
%   That is, the Hamiltonian is H= - sum_k z(k) z(k+1) + B*sum_k x(k), 
%   where x and z denote Pauli spin matrices.

function E0=ising_classical_ground(B)
% See G. Toth, Phys. Rev. A 71, 010301 (R) (2005).
if abs(B)>2,
    E0=-B;
else
    E0=-(1+B^2/4);
end %if

