% keep   Reduced density matrix keeping the given qubits.
%    keep(rho,list) removes the qubits not given in the list.
%    qubits are numbered between 1 and N for an N-qubit state.
%    The expression keep(kron(m2,m1),2) gives m2, while
%    keep(kron(m2,m1),1) gives m1. 
%    (m's are density matrices of pure states.)
%    Thus mkron(m4,m3,m2,m1) would indicate how matrices
%    are numbered for keep.
%    If rho is a state vector, it is converted into
%    a normalized density matrix.
%    If rho is an unnormalized matrix then it is
%    normalized to have trace 1.

function r=keep(rho_in,listneg,varargin)

if isempty(varargin),
   d=2;
else
    if length(varargin)~=1,
        error('Wrong number of input arguments.');
    end %if
    d=varargin{1};
end %if

rho_in=ketbra2(rho_in);

[sx,sy]=size(rho_in);
N=log2(sx)/log2(d);
N=floor(N+0.5);

list=setdiff(N:-1:1,listneg);
%listneg=-sort(-listneg);
Nr=length(listneg);

rho=reorder(rho_in,[listneg list],d);

rho_red=zeros(Nr,Nr);
for k=0:d^Nr-1
   for l=0:d^Nr-1
      rr=0;
      i1=1+k*d^(N-Nr);
      i2=1+l*d^(N-Nr);
      for n=0:d^(N-Nr)-1
         rr=rr+rho(i1+n,i2+n);   
      end %for
      rho_red(k+1,l+1)=rr;
   end %for
end %for
r=rho_red;  
  
  
  