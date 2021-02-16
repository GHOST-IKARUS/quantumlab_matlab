%mkron   Kronecker (tensor) product of several matrices.
%   That is, mkron(X,Y,Z)=kron(kron(X,Y),Z).

function K = mkron(A,B,varargin)
K=kron(A,B);
for n=1:length(varargin)
    K=kron(K,varargin{n});
end %for
