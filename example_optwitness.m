% example_optwitness   Example for the application of example_optwitness
%    We look for a witness for a six-qubit symmetric Dicke state
%    Used in http://arxiv.org/abs/0903.3910.

echo on

% Looking for maximum for Jx^2+Jy^2 for states that are biseparable
% with respect to some bipartition of the qubits

% Number of quibits (2..6)
N=6;

% N-qubit symmetric Dicke state with N/2 excitations
rho=ketbra(dstate(N/2,N));

% White noise
rho_noise=mmstate(N); 

% Projector-based witness
witness_proj=N/(N-1)/2*qeye(N)-rho; 

% Define Pauli spin matrices x,y,z and e
paulixyz;

% Define the collective angular momentum components
Jx=coll(x,N)/2;
Jy=coll(y,N)/2;
Jz=coll(z,N)/2;

% Look for a witness as a linear combination of the followin operators
op_array=[qeye(N) Jx^2 Jx^4 Jx^6 Jy^2 Jy^4 Jy^6 Jz^2 Jz^4 Jz^6]; 

% Determine optimal witness
[W_opt,alpha,coeff]=optwitness(rho,rho_noise,witness_proj,op_array);

echo off

