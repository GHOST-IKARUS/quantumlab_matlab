% Fjinv   Inverse of the function Fj(x) appearing in extreme spin squeezing 
%         of Sorensen and Molmer.
%         See: http://arxiv.org/abs/quant-ph/0011035
%         Only for integer j and x<1/2.
%         Fj_inv(x,j) gives F_j^{-1}(x). The algorithm looks for the 
%         ground state 
%
%         H=Jz^2-lambdax*Jx;
%
%         It searches for the lambdax corresponding to the state 
%         with var(Jz)=xj. Then, it gives back the expectation value of Jx
%         for that state.
%
%      An addtional argument can be given
%      to set the important constants for the numerical algorithm.
%      (If not given, default values ae taken.)
%      In this case, the form is
%
%      Fj_inv(x,j,[lambdaxL,lambdaxH,treshold,lambdatreshold])
%
%      lambdaxL=lower starting value for the Lagrange multiplyier lambdax
%      lambdaxH=higer starting value for the Lagrange multiplyier lambdax
%      treshold=the iteration continues until 
%               (var(Jz)-var(Jz)_input)<treshold j
%      lambdatreshold= search ends with error if
%               (lambdaxH-lambdaxL)<lambdatreshold
%      (In this case the lambdax we need is outside of the interval.)
%      For the default value, see the file.
%
%      See example_Fj_inv.m 

% Program by Geza Toth (C) 2013.

function exJz_over_j=Fj_inv(x,j,varargin)


if length(varargin)==1,
    M=varargin{1};
    lambdaxL=M(1);
    lambdaxH=M(2);
    treshold=M(3);
    lambdatreshold=M(4);
    
elseif length(varargin)==0,
    lambdaxL=0.00000001;
    lambdaxH=20000;
    % Treshold in the accuracy of x
    treshold=0.0000000001;
    % Treshold to test whether lamnbaxL amd lambadxH are almost equal
    % This is the case if the point is outside of the interval
    lambdatreshold=1e-12;
else
    error('Only 2 or 3 arguments are allowed for Fj.')
end %if
if x>1/2
    error('Implemented only for x<=1/2.');
end %if

% The expectation value of jx from the input
% exJx_input=x*j;

% The variance of jz from the input
vaJz_input=x*j;

% k-producabilty (k must be even!)
k=2*j;
%if ~even(k)
if mod(k,2)>0
    error('Implemented only for integer j.');
end %if

% Dimension of matrices describing Jx, Jy and Jz
dd=k+1;
[Jx,Jy,Jz]=su2(dd);

% Initial lambdax
lambdax=(lambdaxH+lambdaxL)/2;

flag=1;

while flag
    
    % Obtain the spin squeezed state as a ground state
    H=Jz^2-lambdax*Jx;
    [v,d]=eig(H);
    DeltaH=0.00001;
    indices=find(real(diag(d))<min(real(diag(d)))+DeltaH);
    if length(indices)>1,
        disp('Degenerate eigenvalues!')
    end %if
    
    % Ground state
    groundstate=v(:,indices(1));
    
    % Expectation value of Jx
    %exJx=ex(Jx,groundstate);
    
    % Variance of Jz
    vaJz=va(Jz,groundstate);
    
    if  vaJz<vaJz_input,
        lambdaxL=lambdax;
        lambdax=(lambdaxH+lambdaxL)/2;
    else
        lambdaxH=lambdax;
        lambdax=(lambdaxH+lambdaxL)/2;
    end %if
    
    difference=vaJz-vaJz_input;
    
    flag=abs(difference)>treshold*j;
    
    if abs(lambdaxL-lambdaxH)<lambdatreshold
       lambdaxH
       lambdaxL
       error('Point outside of the interval. Change lambdaxL and lambdaxH.') 
    end %if
    
end %if

% varJz_over_j=va(Jz,groundstate)/j;
exJz_over_j=ex(Jx,groundstate)/j;


