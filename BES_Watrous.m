% BES_Watrous   Watrous' bound entangled state 
%    The state is defined in
%    https://arxiv.org/abs/quant-ph/0312123 
%    usage: BES_Watrous(epsilon,d) where epsilon and d are the parameters
%    of the state. The density matrix has a dimension d^4 x d^4. 
%    reorder(BES_Watrous(epsilon,d),[4 2 3 1],d) is a state of a 
%    d^2 x d^2 system with interesting distillability properties.

function rho=BES_Watrous(epsilon,d)

F=reordermat([1 2],d);

R=(eye(d^2,d^2)-F)/2;

S=(eye(d^2,d^2)+F)/2;

rho=(d+1+epsilon)/(d-1)*kron(R,R)+kron(S,S);

% test
% difference=abs(trace(rho)-(d+1)*d^3/2-epsilon*d^2*(d-1)/4)

rho=rho/trace(rho);



