% maxeig   Maximum eigenvalue of a matrix 
%    maxeig(M) gives back max(real(eig(M))).
%    Note the function real() in the expression.
%    This takes care of the small imaginary parts
%    appearing in calculations with MATLAB, which could
%    disturb the routine looking for the maximum.

function m=maxeig(M);

m=max(real(eig(M)));

