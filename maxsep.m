% maxsep   Maximum for multi-qudit separable states. 
%   maxsep(op) gives the maximum of op for separable multi-qubit 
%   states. maxsep(op,d) makes it possible to handle qudits
%   with dimension d. Uses simple numerical search.
%   maxsep(op,d,par) makes it possible to set parameters.
%   par is a three-element vector. Defaults value
%   [10000 20000 0.005 ]. First element: Number of
%   random trials in the first phase. Second element:
%   Number of random trials in the second phase.
%   Third element: Constant determining accuracy.
%   Slower than maxsymsep.

function m=maxsep(op,varargin)

% Parameters for the simulation
Delta=0.005;
Nit1=10000;
Nit2=20000;

if isempty(varargin),
   d=2;
elseif length(varargin)==1,
   d=varargin{1};  
elseif length(varargin)==2,
   d=varargin{1};  
   par=varargin{2};
   Nit1=par(1);
   Nit2=par(2);
   Delta=par(3);
else
        error('Wrong number of input arguments.');
end %if

[sx,sy]=size(op);
N=floor(log(sx)/log(d)+0.5);

rmax=-Inf;

fa=zeros(d,N);
famax=fa;
for n=1:Nit1
     %%if mod(n,100)==0,  randn('state',sum(100*clock));  end %if
     f=randn(d,1)+i*randn(d,1);
     fa(:,1)=f;
     for n=1:N-1
        f2=randn(d,1)+i*randn(d,1);
  	 fa(:,n+1)=f2;
        f=kron(f,f2);
     end %for
     r=real(trace(op*f*f')/(f'*f));
     if r>rmax,
         rmax=r;
	 famax=fa;
     end %if
end %for


fa0=famax;
r0=rmax;

% Second phase of the search
for n=1:Nit2
     %%if mod(n,100)==0,  randn('state',sum(100*clock));  end %if
     f=randn(d,1)+i*randn(d,1);
     f=fa0(:,1)+Delta*f;
     fa(:,1)=f;
     for n=1:N-1
        f2=randn(d,1)+i*randn(d,1);
	f2=fa0(:,n+1)+Delta*f2;
	fa(:,n+1)=f2;
        f=kron(f,f2);
     end %for
     r=real(trace(op*f*f')/(f'*f));
     if r>r0,
         r0=r;
         fa0=fa;        
     end %if
end %for
 
m=r0;

