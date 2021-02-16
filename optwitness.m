% optwitness   Optimal entanglement witness for detecting genuine
%              multiqubit entanglement. Uses yalmip and SeDuMi.
%   Form of use: optwitness(rho,rho_noise,witness_proj,op_array), where
%   rho is the state around which we detect entanglement, rho_noise
%   is the noise, witness_proj is the projector witness to the state 
%   around which we detect entanglement and op_array is the array
%   of operators from which we would like to construct the witness.
%   That is, if the three basis operators are B1, B2 and B3, then
%   op_array=[B1,B2,B3]. The form
%   [witness,alpha,coeff]=optwitness(rho,rho_noise,witness_proj,op_array)
%   resurns with the constant alpha for which
%   witness-alpha*witness_proj>=0. coeff contains
%   the coefficients of the basis operators.
%   Invokes the packages SeDuMi and YALMIP.
%   These have to be installed if one wants to use optwitness.
%   See example_optwitness.m.

function [witness,alpha2,coeff2]=optwitness(rho,rho_noise,witness_proj,op_array)

% Extract number of qubits
[sy,sx]=size(rho);
N=log2(sx);
N=floor(N+0.5);

% Number of basis matrices
[sy,sx]=size(op_array);
Nbasismat=sx/sy;

% Variables for minimzation
coeff=sdpvar(Nbasismat,1);
alpha=sdpvar(1,1);

% Compute sum_k coeff_k Basismat_k 
sum_ck_Bk=sdpvar(2^N,2^N);
tr_sum_ck_Bk_rho=sdpvar(1,1);
sum_ck_Bk=coeff(1)*op_array(:,1:2^N);
tr_sum_ck_Bk_rho=coeff(1)*trace(rho*op_array(:,1:2^N));
for k=2:Nbasismat
   sum_ck_Bk=sum_ck_Bk+coeff(k)*op_array(:,(1:2^N)+(k-1)*2^N);
   tr_sum_ck_Bk_rho=tr_sum_ck_Bk_rho+coeff(k)*trace(rho*op_array(:,(1:2^N)+(k-1)*2^N));
end %for

% Constraints
% Important -1>= and not -1== 
F=set(alpha>=0)+set(sum_ck_Bk-alpha*witness_proj>=0)+set(-1>=tr_sum_ck_Bk_rho);
% Solve the semidefinite program
diagnostic=solvesdp(F,trace(sum_ck_Bk*rho_noise));
% Extract the result
witness=double(sum_ck_Bk);

alpha2=double(alpha);
coeff2=double(coeff);

diagnostic

