% isingp   Hamiltonian for the Ising model with periodic BEC
%   isingp(B,N) gives the Hamiltonian for the 
%   ferromagnetic Ising model in transverse field for N qubits, if the 
%   field strength is B and the coefficient 
%   of the nearest neighbor coupling is 1.
%   That is, the Hamiltonian is H= - sum_k z(k) z(k+1) + B*sum_k x(k), 
%   where x and z denote Pauli spin matrices.
%   The coupling is z-z and the direction of the field is x. 
%   For the Hamiltonian peridoc boundary condition is used.
%   If argument n is omitted then the default
%   is taken to be the value of global variable N.

function H=isingp(BField,varargin)

if isempty(varargin),
    global N;
else
    if length(varargin)~=1,
        error('Wrong number of input arguments.')
    end %if
    N=varargin{1};
end %if

x=[0 1;1 0];
y=[0 -i;i 0]; 
z=[1 0; 0 -1];
e=[1 0; 0 1];

% Using routines from the library to make the Ising Hamiltonian
H=-nnchainp(z,z,N)+BField*coll(x,N);
 


