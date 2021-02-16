% example2   Example script for demonstrating the use of QUBIT4MATLAB. 

echo on

% *************************************
% *                                   *
% *          Qubit register           *
% *                                   *
% *************************************

% Define the state (|00>+|11>)|1>/sqrt(2)
psi=ket([0 1 0 0 0 0 0 1])

% Print it out in a pretty form
printv(psi)

% Flip qubits 1 and 2
psi2=reorder(psi,[3 1 2])

% Print it out in a pretty form
printv(psi2)

% Reduced density matrix for psi for qubit 1
rho_red=keep(psi,1)

% The same using remove
rho_red=remove(psi,[3 2])

% Partial transpose according to the third qubit
rho_pt=pt(psi,3)

echo off
