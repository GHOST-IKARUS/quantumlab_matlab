% BES_Horodecki2x4   Horodecki's 2x4 bound entangled state
%    BES_Horodecki2x4(b) gives Horodecki's two-qudit 
%    bound entangled state with the parameter b between 0 and 1.
%    See quant-ph/9801069v1.

function rho=BES_Horodecki2x4(b)

rho=diag([b b b b b b b b]);
rho(5,5)=(1+b)/2;
rho(8,8)=(1+b)/2;
rho(5,8)=sqrt(1-b*b)/2;
rho(8,5)=sqrt(1-b*b)/2;
rho(6,1)=b;
rho(7,2)=b;
rho(8,3)=b;
rho(1,6)=b;
rho(2,7)=b;
rho(3,8)=b;
rho=rho/(7*b+1);



