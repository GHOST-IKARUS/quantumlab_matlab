% pt_nonorm    partial transposition of a density matrix
%              for a qudit register; the density matrix is not
%              normalized
%    pt_nonorm(M,list,d) computes the matrix obtained from M
%    by partially transposing the qudits given in the list.
%    The numbering of the qudits starts with 1, not from 0.
%    Here d is the dimension of the qudits.
%    If d is not provided then d is taken to be 2 (qubits).
%    Unlike for pt, the density matrix is not normalized.

function rhoT=pt_nonorm(rho,list,varargin)

% There is no normalization
% rho=ketbra2(rho);

if nargin==2,
   d=2;
elseif nargin==3,
   d=varargin{1};
else    
   error('Wrong number of input arguments.');
end %if
 
   for n=1:length(list)
      nn=list(n);
      nmax=d^(nn-1);
      nmax2=d^nn;

      % This would be logical:
      %[y,x]=size(rho);
      %rhoT=zeros(y,x);
      % but for using it with yalmip for sdpvar objects, rhoT must
      % be the of the same type. Thus, for intializing a matrix
      % to be of the same size and of the same type as rho,
      % we use the following
      rhoT=rho;

      % Partial transpostion #1
      for k=1:nmax
      for l=1:nmax
          rhoT(k:nmax:end,l:nmax:end)=rho(k:nmax:end,l:nmax:end).';
      end %for
      end %for

      rho=rhoT;

      % Partial transpostion #2
      for k=1:nmax2
      for l=1:nmax2
          rhoT(k:nmax2:end,l:nmax2:end)=rho(k:nmax2:end,l:nmax2:end).';
      end %for
      end %for 
      
      rho=rhoT;
   end %for
   

