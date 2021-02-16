% example_BES_metro4x4   Tests the state BES_metro2018.
%    The state is for a 4x4 system and is given in 
%    Phys. Rev. Lett. 120, 020506 (2018); see
%    https://arxiv.org/abs/1709.03995

clear all
close all
format compact

rho=BES_metro4x4;

sizexy=size(rho)

E=eye(4,4);
D=diag([1 1 -1 -1]);
H=kron(E,D)+kron(D,E);

% The two should be the same
% Quantum Fisher information obtained numerically
fisher_rho_ppt=fisher(rho,H)

% Quantum Fisher information obtained from the article
fisher_rho_ppt_theory=16*sqrt(2)/(1+sqrt(2))

% Variance
var_H_rho_ppt_times4=va(H,rho)*4

% The minum eigenvalue of the partial transpose should be non-negative
mineig_rho_pt=mineig(pt(rho,1,4))




