% mestate   Maximally entangled state
%    mestate(d) gives the state vector of a maximally 
%    entangled state of two qudits of dimension d.

function s=mestate(d)

s=zeros(d*d,1);
for n=0:d-1
    s(n+d*n+1)=1;
end %for
s=s/sqrt(s'*s);