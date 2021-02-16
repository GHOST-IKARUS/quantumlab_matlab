% nm   Normalization
%    nm(v) converts vector v into a column vector and
%    normalizes it. If v is not a vector, then
%    it normalizes it as a density matrix setting
%    the trace to 1.

function w=nm(v)
   [y,x]=size(v);
   if y==1,
     w=v.';
     w=w/sqrt(w'*w);    
   elseif x==1,
     w=v;
     w=w/sqrt(w'*w);
   else
     w=v/trace(v);
   end %if
