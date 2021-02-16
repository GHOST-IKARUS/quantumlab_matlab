%   maxb(op) gives maximum value for
%   operators for biseparable states; 
%   considers states with all possible
%   bipartitioning. Uses simple numerical search.
%   maxb(op,par) makes it possible to set parameters.
%   par is a three-element vector. Defaults value
%   [10000 20000 0.005 ]. First element: Number of
%   random trials in the first phase. Second element:
%   Number of random trials in the second phase.
%   Third element: Constant determining accuracy.
%   Calls maxbisep.m.

function mb=maxb(op1,varargin)

if isempty(varargin),
    % Parameters for the simulation
    par=[10000 20000 0.005];
elseif length(varargin)==1,
   par=varargin{1};
else
   error('Wrong number of input arguments.');
end %if

% Qubits
d=2;

[sx,sy]=size(op1);
N=log2(sx)/log2(d);
N=floor(N+0.5);

% Loop for all bipartitions
mb=0;
for n=1:2^N-2
    list=[];
    b=dec2bin(n+2^N);
    for m=N:-1:1
       if b(m+1)=='1',
           list=[list,m];
       end %if
    end %for 
    if length(list)<=floor(N/2+0.01),
       mb=max(mb,max(eig(maxbisep(op1,list,par))));
   end %if
end %for