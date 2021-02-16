%nnchain   Defines a Hamiltonian with a(k)b(k+1) nearest-neighbor
%          interaction
%   nnchain (a,b,n) defines a a(k)b(k+1) type Hamiltonian
%   with non-periodic boundary condition, n is the number of qudits.
%   If argument n is omitted than the default is taking to be
%   the value of global variable N.
%   Works even for chains of qudits with a dimension larger than two.

function h=nnchain(a,b,varargin)
if isempty(varargin),
    global N;
else
    if length(varargin)~=1,
        error('Wrong number of input arguments.')
    end %if
    N=varargin{1};
end %if
% Obtain dimension
[d,junk]=size(a);
h=zeros(d^N);
op=kron(a,b);
for n=1:N-1
    h=h+kron(kron(eye(d^(n-1)),op),eye(d^(N-n-1)));
end %for
