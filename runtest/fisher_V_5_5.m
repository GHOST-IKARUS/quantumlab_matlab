% fisher   Quantum Fisher information. Usage: fisher(rho,A), where
%          rho is a density matrix and A is an operator.
%          Alternatively, it can also be used as fisher(rho,A,B).
%          Note that fisher(rho,A,A)=fisher(rho,A).
%          The form fisher(rho,A,B,threshold) makes it possible
%          to define the threshold below which an eigenvalue is
%          consdiered zero. The defaults value is 1e-20.

function f=fisher(rho,J,varargin)

% Make rho surely Hermitian
% A little bit of non-Hermitianicity can lead to
% non-orthogonal eigenvalues for eig
rho=(rho+rho')/2;

%Threshold=1e-10;
Threshold=1e-20;

if length(varargin)==0

   [v,d]=eig(rho);
   [sy,sx]=size(rho);
   f=0;
   for n=1:sx
       for m=1:n-1
           lambdan=d(n,n);
           lambdam=d(m,m);
           if abs(lambdam+lambdan)>Threshold;
               f=f+(lambdam-lambdan)^2/(lambdam+lambdan)*abs(braket(v(:,n),J,v(:,m)))^2;
           end %if      
       end %for
   end %for

   f=f*4; % Because of missing factor of 2, and because we count only half of the terms  
   
else
    
   if length(varargin)==2
      Threshold=varargin{2};
   end %if
      
   J2=varargin{1};
    
   [v,d]=eig(rho);
   [sy,sx]=size(rho);
   f=0;
   for n=1:sx
       for m=1:n-1
           lambdan=d(n,n);
           lambdam=d(m,m);
           if abs(lambdam+lambdan)>Threshold;
               f=f+(lambdam-lambdan)^2/(lambdam+lambdan)*...
               (braket(v(:,n),J,v(:,m)))*(braket(v(:,m),J2,v(:,n)));
           end %if      
       end %for
   end %for

   f=f*4; % Because of missing factor of 2, and because we count only half of the terms
   
end %if
    
