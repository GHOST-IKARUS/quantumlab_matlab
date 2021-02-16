%dstate   Defines a symmetric Dicke state.
%   dstate(e,n) gives the state vector for an n-qubit symmetric
%   Dicke state with e excitations. If argument n is omitted 
%   then the default is taken to be the value of global variable N.

function w=dstate(e,varargin)
if isempty(varargin),
    global N;
else
    if length(varargin)~=1,
        error('Wrong number of input arguments.')
    end %if
    N=varargin{1};
end %if
w=zeros(2^N,1);
for n=1:2^N
   % Conversion to binary 
   b=dec2bin(n-1+2^N);
   b=b(2:end);  
   a=[];
   for k=1:length(b)
      a=[a;str2num(b(k))];
   end %for
   if sum(a)==e,
      w(n)=1;
   end %if
end %for
w=w/sqrt(w'*w);
