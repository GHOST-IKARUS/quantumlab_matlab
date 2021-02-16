% example_maxppt   Example for the application of maxppt
%    This example computes the maximum for the operator Jx^2+Jy^2 
%    for states that are PPT with respect to some bipartition of the qubits.
%    Similar otimisation is carried out in http://arxiv.org/abs/0903.3910.

echo on

% Looking for maximum for Jx^2+Jy^2 for states that are biseparable
% with respect to some bipartition of the qubits

% Number of qubits 
N=4;

% Define Pauli spin matrices x,y,z and e
paulixyz;

% Define the collective angular momentum components
Jx=coll(x,N)/2;
Jy=coll(y,N)/2;
Jz=coll(z,N)/2;

% Define the operator for which we look for the maximum
M=Jx^2+Jy^2;

% Maximum for the 1:234 partition
max1=maxppt(M,1)

% Maximum for the 12:34 partition
max12=maxppt(M,1:2)

echo off

