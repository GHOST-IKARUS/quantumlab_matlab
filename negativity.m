% negativity   Negativity for a qudit register.
%    negativity(M,list,d) computes the matrix obtained from M
%    by partially transposing the qudits given in the list,
%    and returns with the sum of the absolute value of the 
%    negative eigenvalues of this matrix. 
%    Note that M is normalized to be a density matrix,
%    before computing the eigenvalues.
%    That is, it computes the negativity
%    corresponding to a given partitioning.
%    The numbering of the qudits starts with 1.
%    Here d is the dimension of the qudits.
%    If d is not provided then it is taken to be 2 (qubits).

function neg=negativity(rho,list,varargin);

if isempty(varargin),
   d=2;
else
    if length(varargin)~=1,
        error('Wrong number of input arguments.');
    end %if
    d=varargin{1};
end %if

% Convert state vector to density matrix if necessary
rho=ketbra2(rho);

[v,d]=eig(pt(rho,list,d));

d=diag(d);

% real is needed because MATLAB tends to produce
% small imaginary parts in some computation
d=real(d);

neg=-sum(d.*(d<0));


