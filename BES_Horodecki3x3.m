% BES_Horodecki3x3   Horodecki's 3x3 bound entangled state
%    BES_Horodecki3x3(a) gives Horodecki's two-qutrit 
%    bound entangled state with the parameter a between 0 and 1.
%    See quant-ph/9801069v1.

function rho=BES_Horodecki3x3(a)

rho=diag([a a a a a a a a a]);
rho(7,7)=(1+a)/2;
rho(9,9)=(1+a)/2;
rho(9,7)=sqrt(1-a*a)/2;
rho(7,9)=sqrt(1-a*a)/2;
rho(5,1)=a;
rho(9,1)=a;
rho(5,9)=a;
rho(1,5)=a;
rho(1,9)=a;
rho(9,5)=a;
rho=rho/(8*a+1);



