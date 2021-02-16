%gstate_stabilizer    Gets the generators for the stabilizer of a graph state.
%   gstate_stabilizer(Gamma) gives the generator for the stabilizer 
%   for a graph state corresponding to a connectivity matrix Gamma.

function stab=gstate_stabilizer(Gamma)

x=[0 1;1 0];
z=[1 0;0 -1];
y=i*x*z;

[sy,sx]=size(Gamma);
N=sy;

stab=cell(N,1);
for n=1:N
   op=1;
   for m=1:n-1
      if Gamma(n,m)>0,
         op2=z;
      else
         op2=eye(2);
      end %if
      op=kron(op,op2);    
   end %for
   op=kron(op,x);
   for m=n+1:N
      if Gamma(n,m)>0,
         op2=z;
      else
         op2=eye(2);
      end %if
      op=kron(op,op2);    
   end %for
   stab{n}=op;
end %for


