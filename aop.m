% aop   Gives the annihilation operator as a matrix. Usage: aop(d),
%       where d is the size of the matrix.

function m=aop(d);

m=zeros(d,d);
for n=1:d-1
    m(n,n+1)=sqrt(n);    
end %for