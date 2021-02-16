% mrealign   Generalized realignment of a multipartite operator
%    mrealign(M,list) realigns matrix M (that can be a density matrix)
%    in a multipartite system where the qudits have dimension 2.
%    mrealign(M,list,d) makes it possible to have qudits
%    for d>2. 
%    Example for 2 qubits:
%       mrealign(M,[4 3 2 1]) - Noting changes
%       mrealign(M,[2 3 4 1]) - Partial transpose
%       mrealign(M,[4 2 3 1]) - Realignment

function rhoR=mrealign(rho,list,varargin)

% rho is should not be normalized
% thus we do not use: "rho=ketbra2(rho);"

if isempty(varargin),
   d=2;
else
    if length(varargin)~=1,
        error('Wrong number of input arguments.');
    end %if
    d=varargin{1};
end %if

[sx,sy]=size(rho);
N=floor(log(sx)/log(d)+0.5);

% Create a multidimensional matrix
r2=reshape(rho,kron(d,ones(1,2*N)));

% Permute the indices
% Nothing happens corresponds to [2N 2N-1 .... 3 2 1]
list2=2*N+1-list;
r3=permute(r2,list2);
  
rhoR=reshape(r3,d^N,d^N);
  