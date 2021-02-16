% cohstate   cohstate(alpha,M) gives the state vector of M elements for 
%            a coherent state with given alpha 

function v=cohstate(alpha,M)

v=zeros(M,1);
v((0)+1)=1;
for n=0:M-2
    v((n+1)+1)=alpha/sqrt(n+1)*v((n)+1);   
end %for
v=v/sqrt(v'*v);
