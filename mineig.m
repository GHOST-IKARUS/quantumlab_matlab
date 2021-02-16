% mineig   Minimum eigenvalue of a matrix 
%    mineig(M) gives back min(real(eig(M))).
%    Note the function real() in the expression.
%    This takes care of the small imaginary parts
%    appearing in calculations with MATLAB, which could
%    disturb the routine looking for the minimum.

function m=mineig(M);

m=min(real(eig(M)));

