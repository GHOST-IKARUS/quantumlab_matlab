% ketbra2   Creating a density matrix from an unnormalized state vector.
%    Creates a density matrix from a state vector. If it was
%    already a density matrix then it just normalizes it.

function k=ketbra2(f);
[sy,sx]=size(f);
if sx==1,
   k=f*f';
elseif sy==1,
   k=f'*f;
else
   k=f;
end %if
if trace(k)~=0,
   k=k/trace(k);
end %if
