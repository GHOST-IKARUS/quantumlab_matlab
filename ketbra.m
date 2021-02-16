% ketbra   Creating a density matrix from an unnormalized state vector.
%    Creates a density matrix from a state vector. If it was
%    not a column or a row vector then it sends an error message.

function k=ketbra(f);
[sy,sx]=size(f);
if sx==1,
   k=f*f';
elseif sy==1,
   k=f'*f;
else
   error('Row or column vector expected.')
end %if
if trace(k)~=0,
   k=k/trace(k);
end %if
