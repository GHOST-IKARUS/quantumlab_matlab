%gstate    Defines a graph state.
%   gstate(Gamma) gives the state vector for a graph state
%   corresponding to a connectivity matrix Gamma.
%   The form [g,stab]=gstate(Gamma) gives the state
%   vector for the graph state in vector g and the 
%   generator of the stabilizer in the cell array stab.

function [g,stab]=gstate(Gamma)

x=[0 1;1 0];
z=[1 0;0 -1];
y=i*x*z;

% Get the generators of the stabilizer
stab=gstate_stabilizer(Gamma);

W=0;
for n=1:length(stab)
   W=W-stab{n};
end %for n

[v,d]=eig(W);
[junk,index]=min(diag(d));
g=v(:,index(1));
g=g/sqrt(g'*g);

