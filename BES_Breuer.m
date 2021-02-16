% BES_Breuer   Breuer's two-qudit bound entangled state 
%              for dimension d>=4, even d
%    BES_Breuer(lambda,d) gives Breuer's bound entangled state
%    for two qudits of dimension d, with the lambda paramater
%    describing the weight of the singlet component
%    For the state see H.-P. Breuer, Phys. Rev. Lett. 97, 080501 (2006).
%    Journal article: http://link.aps.org/doi/10.1103/PhysRevLett.97.080501
%    Preprint:        http://arxiv.org/abs/quant-ph/0605036

function rho=BES_Breuer(lambda,d)

% SU(2) generators for dxd matrices
N=d-1;
a1da1=diag([0:N]);
a2da2=diag([N-(0:N)]);
a1da2=zeros(N+1,N+1);
for k=0:N-1
    a1da2(k+1+1,k+1)=sqrt((k+1)*(N-k));   
end %for
a2da1=a1da2';
% Schwinger's construction
jz=(a1da1-a2da2)/2;
jx=(a1da2+a2da1)/2;
jy=-i*(a1da2-a2da1)/2;
% Collective operators
ee=eye(d);
Jx=kron(jx,ee)+kron(ee,jx);
Jy=kron(jy,ee)+kron(ee,jy);
Jz=kron(jz,ee)+kron(ee,jz);

% A singlet is the ground state of the following Hamiltonian
H=Jx^2+Jy^2+Jz^2;
[v,D]=eig(H);
[smallestE,index]=min(diag(real(D)));
phis=ket(v(:,index));

% Mix it with the normalized projector to the symmetric
% subspace
rho=lambda*ketbra(phis)+(1-lambda)*nm(proj_sym(2,d));


