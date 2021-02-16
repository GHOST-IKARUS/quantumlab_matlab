% optspinsq(rho)   Optimal spin squeezing inequalities
%    optspinsq(rho) gives back a negative value if the
%    multi-qubit state rho is detected as entangled
%    by the optimal spin squeezing inequalities.
%    (See Eq. (2.b,c,d) http://arxiv.org/abs/quant-ph/0702219.
%     Eq. (2.a) is not included since it is satisfied by any quantum state.) 
%    The form [fmin,f123]=optspinsq(rho) gives back
%    in f123 a three element array. Each element of the
%    array gives -1 times the violation of the corresponding spin
%    squeezing inequality. fmin is the minimum of the three values.
%    If one of them is negative then the state is detected as entangled.
%    Beside the inequalities themselves, a method is also
%    implemented that looks for the optimal choice of x, y, and z
%    coordinates (See fourth page of paper above).

function [fmin,f123]=optspinsq(rho)

[sy,sx]=size(rho);
N=log2(sx);

% Collective observables
paulixyz;
Jx=coll(x,N); % Here I did not use the 1/2!
Jy=coll(y,N);
Jz=coll(z,N);

% Correlation matrix
Cxx=trace(rho*Jx*Jx);
Cyy=trace(rho*Jy*Jy);
Czz=trace(rho*Jz*Jz);
Cxy=trace(rho*(Jx*Jy+Jy*Jx)/2);
Cxz=trace(rho*(Jx*Jz+Jz*Jx)/2);
Cyz=trace(rho*(Jy*Jz+Jz*Jy)/2);
C=[Cxx Cxy Cxz;Cxy Cyy Cyz; Cxz Cyz Czz];

% Covariance matrix
v=[trace(Jx*rho) trace(Jy*rho) trace(Jz*rho)];
g=C-v.'*v;

X=(N-1)*g+C;
% /4 is needed since the definition of Jx (Jx=coll(x,N)) did not include a
% factor of 1/2
f123=[trace(g)-2*N,mineig(X)-trace(C)+2*N,(N-1)*trace(g)-N*(N-2)-maxeig(X)]/4;
fmin=min(real(f123));



