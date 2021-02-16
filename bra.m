% bra  After element-wise conjugation, transforms a vector into a normalized row vector.

function w=bra(v)
   [y,x]=size(v);
   if x>1,
     w=conj(v);
   else
     w=v';
   end %if
   % normalization
   w=w/sqrt(w*w');
   
   