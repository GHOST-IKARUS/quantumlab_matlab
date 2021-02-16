% Fj   The function Fj(x) appearing in extreme spin squeezing 
%      of Sorensen and Molmer.
%      See: http://arxiv.org/abs/quant-ph/0011035
%      Only for integer j and x<1.
%      Fj(x,j) gives F_j(x). The algorithm looks for the ground state 
%
%         H=Jz^2-lambdax*Jx;
%
%      It searches for the lambdax corresponding to the state 
%      with ex(Jx)=xj. Then, it gives back var(Jz) for that state.
%
%      An addtional argument can be given
%      to set the important constants for the numerical algorithm.
%      (If not given, default values ae taken.)
%      In this case, the form is
%
%      Fj(x,j,[lambdaxL,lambdaxH,treshold,lambdatreshold])
%
%      lambdaxL=lower starting value for the Lagrange multiplyier lambdax
%      lambdaxH=higer starting value for the Lagrange multiplyier lambdax
%      treshold=the iteration continues until 
%               (ex(Jx)-ex(Jx)_input)<treshold j
%      lambdatreshold= search ends with error if
%               (lambdaxH-lambdaxL)<lambdatreshold
%      (In this case the lambdax we need is outside of the interval.)
%      For the default value, see the file.
%
%      See example_Fj.m 

% Program by Geza Toth (C) 2013.

function varJz_over_j=Fj(x,j,varargin)


if length(varargin)==1,
    M=varargin{1};
    lambdaxL=M(1);
    lambdaxH=M(2);
    treshold=M(3);
    lambdatreshold=M(4);
    
elseif length(varargin)==0,
    lambdaxL=0.00000001;
    lambdaxH=200000;
    % Treshold in the accuracy of x
    treshold=0.00000000001;
    % Treshold to test whether lamnbaxL amd lambadxH are almost equal
    % This is the case if the point is outside of the interval
    lambdatreshold=1e-12;
else
    error('Only 2 or 3 arguments are allowed for Fj.')
end %if
if x>1
    error('Implemented only for x<=1.');
end %if

% The expectation value of jx from the input
exJx_input=x*j;

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

% Treshold to detect degeneracy
TresholdH=0.0001;

while flag
    
    % Obtain the spin squeezed state as a ground state
    H=Jz^2-lambdax*Jx;
    [v,d]=eig(H);
    indices=find(real(diag(d))<min(real(diag(d)))+TresholdH);
    if length(indices)>1,
        disp('Degenerate eigenvalues!')
    end %if
    
    % Ground state
    groundstate=v(:,indices(1));
    
    % Expectation value of Jx
    exJx=ex(Jx,groundstate);
    
    if  exJx<exJx_input,
        lambdaxL=lambdax;
        lambdax=(lambdaxH+lambdaxL)/2;
    else
        lambdaxH=lambdax;
        lambdax=(lambdaxH+lambdaxL)/2;
    end %if
    
    difference=exJx-exJx_input;
    
    flag=abs(difference)>treshold*j;
    
    if abs(lambdaxL-lambdaxH)<lambdatreshold
       lambdaxH
       lambdaxL
       error('Point outside of the interval. Change lambdaxL and lambdaxH.') 
    end %if
    
end %if

varJz_over_j=va(Jz,groundstate)/j;


