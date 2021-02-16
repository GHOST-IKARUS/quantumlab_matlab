% orthogobs   Orthogonal observables for a qudit
%   orthogobs(d) gives the array of d^2 local orthogonal observables
%   for a d-state system. obs=orthogobs(d) gives back a d x d x d^2
%   dimensional matrix. The kth observable can be accessed as obs(:,:,k).
%   These observables fulfill trace(A_k*A_l)=0 for k>l and
%   trace(A^2)=1. For the definition of such observables
%   see http://www.arxiv.org/abs/quant-ph/0412220v2.

function obs=orthogobs(d)

obs=zeros(d,d,d^2);
index=1;

for n=1:d    
    v=zeros(d,1);
    v(n)=1;
    obs(:,:,index)=v*v.';
    index=index+1;
end %for

for n=1:d
    for m=1:n-1
        v=zeros(d,1);
        v(n)=1;
        u=zeros(d,1);
        u(m)=1;
        obs(:,:,index)=(u*v.'+v*u.')/sqrt(2);
        index=index+1;
    end %for
end %for

for n=1:d
    for m=1:n-1
        v=zeros(d,1);
        v(n)=1;
        u=zeros(d,1);
        u(m)=1;
        obs(:,:,index)=(u*v.'-v*u.')/i/sqrt(2);
        index=index+1;
    end %for
end %for


