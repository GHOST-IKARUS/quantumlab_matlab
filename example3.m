% example3   Example script for demonstrating the use of QUBIT4MATLAB. 

echo on

% *************************************
% *                                   *
% *    Operators for spin chains      *
% *                                   *
% *************************************

% Define Pauli spin matrices x,y,z and e
paulixyz

% Define Heisenberg interaction for two qubits
H_H=kron(x,x)+kron(y,y)+kron(z,z)

% Print out the decomposition of H_H
decompose(H_H)

% Print out the decomposition of H_H in LaTeX format
decompose(H_H,1)

% Define a spin chain Hamiltonian 
% with Heisenberg interaction
% periodic boundary condition, 8 qubits
% Compute ground state energy.
H_Hp=nnchainp(x,x,8)+nnchainp(y,y,8)+nnchainp(z,z,8);
ground_state_energy=min(real(eig(H_Hp)))

% Define a similar Hamiltonian, but not with a periodic
% boundary condition
H_H=nnchain(x,x,8)+nnchain(y,y,8)+nnchain(z,z,8);
ground_state_energy=min(real(eig(H_H)))

% Define an operator acting on the 3rd and 4th qubit in a 
% 5-qubit chain
Op=twoquditop(kron(x,x)+kron(y,y)+kron(z,z),3,4,5);

% Print out the decompostion of Op
decompose(Op)

echo off