% Convex roof of the linear entropy - example
% compare with CORONA package
% The CORONA package can be found here:
%  http://de.mathworks.com/matlabcentral/fileexchange/47823-corona-convex-roof-numerical-analysis?requestedDomain=www.mathworks.com


clear all
clear classes

rho=BES_Horodecki3x3(0.28);

rank_rho=rank(rho)
e_lin_QUBIT4MATLAB=elin(rho)

dqudit=3;
h = Elin(rho, dqudit, dqudit, 'complex');
elin_CORONA = h.lb_convex_roof()



