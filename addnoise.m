% addnoise   Adds white noise to a density matrix.
%   addnoise(rho,p) computes the matrix 
%   rho'=(1-p)*rho+p*eye(M)/M where rho is an MxM matrix
%   If rho is a state vector then 
%   it is converted into a normalized density matrix.

function rho_noisy=addnoise(rho,p)

% Convert state vector to density matrix if necessary
rho=ketbra2(rho);

% Obtain the size of the density matrix
[sx,sy]=size(rho);
rho_noisy=(1-p)*rho+p*eye(sx)/sx;