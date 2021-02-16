% example_BES_private.   Tests the state BES_private.
%    The state is given in https://arxiv.org/abs/2002.12409
%    and in https://arxiv.org/abs/1309.7992
%    States are handled both in the ABA'B' and the AA'BB' cases.

clear all
close all
format compact

% Set the size of the state, ie, D=4,6,8,10,12,etc
D=4;

% Try the state using the ABA'B' partitioning of the syste, which has been used in the papers.

disp(['ABA' char(39) 'B' char(39) ' order'])

d=D/2;
rho=BES_private(D);

paulixyz;
H=mkron(z,e,eye(d,d),eye(d,d))+mkron(e,z,eye(d,d),eye(d,d));

% Quantum Fisher information obtained numerically
fisher_rho_ppt=fisher(rho,H)

% Quantum Fisher information obtained from the article
% Should be the same as the nuemrical value
fisher_rho_ppt_theory=16*sqrt(d)/(1+sqrt(d))

% The minum eigenvalue of the partial transpose should be non-negative
mineig_rho_pt=mineig(pt(rho,[4 2],[2 2 d d]))

%
% Try the state using the AA'BB' partitioning of the system
%

disp(' ')
disp(['AA' char(39) 'BB' char(39) ' order'])

rho2=BES_private(D,1);
H2=mkron(z,eye(d,d),e,eye(d,d))+mkron(e,eye(d,d),z, eye(d,d));

% Quantum Fisher information obtained numerically
fisher_rho_ppt=fisher(rho2,H2)

% Quantum Fisher information obtained from the article
% Should be the same as the nuemrical value
fisher_rho_ppt_theory=16*sqrt(d)/(1+sqrt(d))

% The minum eigenvalue of the partial transpose should be non-negative
mineig_rho_pt=mineig(pt(rho2,1,2*d))








