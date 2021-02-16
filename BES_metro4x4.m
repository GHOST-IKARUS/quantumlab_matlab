% BES_metro4x4   Bound entangled state for robust metrology
%    The state is for a 4x4 system and is given in 
%    Phys. Rev. Lett. 120, 020506 (2018); see
%    https://arxiv.org/abs/1709.03995

function rho=BES_metro4x4

v0=[1 0 0 0]';
v1=[0 1 0 0]';
v2=[0 0 1 0]';
v3=[0 0 0 1]';

Psi1=(kron(v0,v1)+kron(v2,v3))/sqrt(2);
Psi2=(kron(v1,v0)+kron(v3,v2))/sqrt(2);
Psi3=(kron(v1,v1)+kron(v2,v2))/sqrt(2);
Psi4=(kron(v0,v0)-kron(v3,v3))/sqrt(2);
Psi5=(+kron(v0,v3)+kron(v1,v2))/2+kron(v2,v1)/sqrt(2);
Psi6=(-kron(v0,v3)+kron(v1,v2))/2+kron(v3,v0)/sqrt(2);
q=(sqrt(2)-1)/2;
p=(1-2*q)/4;
rho=p*(Psi1*Psi1'+Psi2*Psi2'+Psi3*Psi3'+Psi4*Psi4')...
   +q*(Psi5*Psi5'+Psi6*Psi6');

