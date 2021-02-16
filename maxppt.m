% maxppt   Maximum for multi-qudit states with a positive partial transpose.
%   maxppt (op,list,d) gives maximum value for an
%   operator op for states with a positive partial transpose
%   for a given bipartitioning. list contains the list of 
%   qubits that are in one group. Qubits are numbered from 1 to N.
%   d is the dimension of qudits. If d is omitted, it is
%   assumed that d=2 (qubits).
%   maxppt uses semidefinite programming, 
%   and invokes the packages SeDuMi and YALMIP.
%   These has to be installed if one wants to use maxppt.
%   Thus, these two packages must be installed for using maxppt.
%   See example_maxppt.m.

function m=maxppt(op,list,varargin)

if nargin==2,
   d=2;
elseif nargin==3,
   d=varargin{1};
else    
   error('Wrong number of input arguments.');
end %if

[sx,sy]=size(op);
Size=sy;
N=floor(log(Size)/log(d)+0.5);

op1=reorder(op,[list setdiff(N:-1:1,list)]);

% Semidefinite program with yalmip
% Define rho, the density matrix
rho=sdpvar(Size,Size,'hermitian','complex');

% Computing the partial transpose
% Partial transpose for the (12..D):(D+1..N) partition
D=d^(N-length(list));
R=1:Size/D;
W=Size/D;
rhoT=[rho(R+0*W,R+0*W).'];
for n=1:D-1
   rhoT=[rhoT,rho(R+0*W,R+n*W).'];
end %for
for m=1:D-1
    A=[rho(R+m*W,R+0*W).'];
    for n=1:D-1
           A=[A,rho(R+m*W,R+n*W).'];
    end %for
    rhoT=[rhoT;A];
end %for

% Solve the semidefinite program
% Important: write rhoT+rhoT'>=0, instead of
% rhoT>=0, otherwise it could count as an elementwise
% condition, if not perfectly symmetric
diagnostic=solvesdp(set(rho>=0)+set(rhoT+rhoT'>=0)+set(trace(rho)==1),-trace(op1*rho));
% Extract the result
m=double(trace(op1*rho));



