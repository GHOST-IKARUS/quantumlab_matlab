% maxbisep   Maximum for biseparable multi-qubit states.
%   maxbisep(op,list) gives maximum value for an
%   operator for biseparable states; 
%   list determines the biparitioning.
%   Uses simple numerical search.
%   maxbisep(op,list,par) makes it possible to set parameters.
%   par is a three-element vector. Defaults value
%   [10000 20000 0.005 ]. First element: Number of
%   random trials in the first phase. Second element:
%   Number of random trials in the second phase.
%   Third element: Constant determining accuracy.

function m=maxbisep(op1,list,varargin)

if isempty(varargin),
    % Parameters for the simulation
    Delta=0.005;
    Nit1=10000;
    Nit2=20000;    
elseif length(varargin)==1,
   par=varargin{1};
   Nit1=par(1);
   Nit2=par(2);
   Delta=par(3);
else
   error('Wrong number of input arguments.');
end %if


% Qubits
d=2;

[sx,sy]=size(op1);
N=log2(sx)/log2(d);
N=floor(N+0.5);

listneg=setdiff(N:-1:1,list);
listneg=-sort(-listneg);
Nr=length(listneg);
if isequal([listneg list],N:-1:1),
  op=op1;
else
  op=reorder(op1,[listneg list],d);
end %if

% Dimensions
d1=d^length(listneg);
d2=d^length(list);

rmax=-Inf;

for n=1:Nit1
     %%if mod(n,100)==0,  randn('state',sum(100*clock));  end %if
     f1=randn(d1,1)+i*randn(d1,1);
     f2=randn(d2,1)+i*randn(d2,1);
     f=kron(f1,f2);
     r=real((f'*op*f)/(f'*f));
     if r>rmax,
         rmax=r;
         f1max=f1;
         f2max=f2;        
     end %if
 end %for
 
 f10=f1max;
 f20=f2max;
 r0=rmax;
 % inital_guess=r0;
 
 % Second phase of the search
 for n=1:Nit2
     %%if mod(n,100)==0,  randn('state',sum(100*clock));  end %if
     f1=randn(d1,1)+i*randn(d1,1);
     f2=randn(d2,1)+i*randn(d2,1);
     f1=f10+Delta*f1;
     f2=f20+Delta*f2;
     f=kron(f1,f2);
     r=real((f'*op*f)/(f'*f));
     if r>r0,
         r0=r;
         f10=f1;
         f20=f2;        
     end %if
end %for
 
m=r0;
% save
% s=kron(f10,f20);

