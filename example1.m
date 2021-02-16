% example1   Example script for demonstrating the use of QUBIT4MATLAB. 

echo on

% *************************************
% *                                   *
% *  A single two-state state system  *
% *                                   *
% *************************************

% Define a state vector corresponding to |0>
phi0=ket([1 0])

% Define a state vector corresponding to |1>
phi1=ket([0 1])

% Define a state vector corresponding to (|0>+|1>)/sqrt(2)
phi01=ket([1 1])

% Density matrix corresponding to a state vector
rho=ketbra(phi01)

echo off
