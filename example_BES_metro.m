% example_BES_metro   Tests the state BES_metro.
%    The state is given in https://arxiv.org/abs/2002.12409
%    States are handled both in the ABA'B' and the AA'BB' cases.

clear all
close all
format compact


% Dimension of the example, can be D=6 or power of 2 such as 4,8,16,32
D=6

% Try the state using the ABA'B' partitioning of the syste, which has been used in the papers.

disp(['ABA' char(39) 'B' char(39) ' order'])

d=D/2;
rho=BES_metro(D);

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


d=D/2;
rho=BES_metro(D,1);

paulixyz;
H=mkron(z,eye(d,d),e,eye(d,d))+mkron(e,eye(d,d),z,eye(d,d));


% Quantum Fisher information obtained numerically
fisher_rho_ppt=fisher(rho,H)

% Quantum Fisher information obtained from the article
% Should be the same as the nuemrical value
fisher_rho_ppt_theory=16*sqrt(d)/(1+sqrt(d))

% The minum eigenvalue of the partial transpose should be non-negative
mineig_rho_pt=mineig(pt(rho,1,2*d))








