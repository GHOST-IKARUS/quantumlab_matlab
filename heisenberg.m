% heisenberg   Hamiltonian for the Heisenberg model with non-periodic BEC
%   heisenberg(N) gives the Hamiltonian for the 
%   anti-ferromagnetic Heisenberg chain for N qubits, if the coefficient 
%   of the nearest neighbor coupling is 1.
%   For the Hamiltonian non-periodic boundary condition is used.
%   If argument N is omitted then the default
%   is taken to be the value of global variable N.

function H=heisenberg(varargin)

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
   H=nnchain(z,z,N)+nnchain(x,x,N)+nnchain(y,y,N);
 


