% Convex roof of the linear entropy - example

clear all

tic

rho=BES_UPB3x3;

rank_rho=rank(rho)
% The result is 0.065191337232140 
E_lin_UPB3x3=elin(rho)


rho=BES_Horodecki3x3(0.28);

rank_rho=rank(rho)
% The results is 0.001666593232004 
E_lin_Horodecki3x3=elin(rho)

toc


