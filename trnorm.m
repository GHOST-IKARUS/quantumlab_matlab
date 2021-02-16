% trnorm   Tracenorm of a matrix

function t=trnorm(A)
   t=trace((A*A')^0.5);
   % Get rid of small imaginary part due to the mysteries of MATLAB
   t=real(t);
