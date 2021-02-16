% ising_free - Free energy of the Ising model in thermal equilibrium
%   ising_free(B,T) gives the free energy of 
%   Ising model in transverse field in the
%   thermodynamic limit, if the temperature is T, the  
%   field strength is B and the coefficient 
%   of the nearest neighbor coupling is 1.

function F=ising_free(BField,T)

J=4;
Gamma=2*BField;
lambda=J/2/Gamma;

k=1;
dq=0.00001;
q=0:dq:pi;

wq=sqrt(1+2*lambda*cos(q)+lambda^2);

delta=log(cosh(0.5/k/T*Gamma*wq));
F=-k*T*(log(2)+1/pi*sum(delta)*dq);


