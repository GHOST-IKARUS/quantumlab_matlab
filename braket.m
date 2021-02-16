% braket  Dirac's bra-ket
%    braket(phi1,phi2) denotes the scalar product of phi1 and phi2.
%    braket(phi1,op,phi2) is the same as bra(phi1)*op*ket(phi2).

function b=braket(phi1,varargin)
   if length(varargin)==1,
       b=bra(phi1)*ket(varargin{1});
   else
       b=bra(phi1)*varargin{1}*ket(varargin{2});
   end %if