% ket  Transforms a vector into normalized column vector.

function w=ket(v)
   [y,x]=size(v);
   if x>1,
     w=v.';
   else
     w=v;
   end %if
   % normalisation
   w=w/sqrt(w'*w);
   
   
