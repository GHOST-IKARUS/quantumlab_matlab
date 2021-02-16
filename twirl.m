% twirl   Twirling
%   twirl(rho) twirls the multi-qubit density matrix rho.
%   [rho,difference]=twirl(rho) gives also the 
%   norm of the difference between the original and the 
%   twirled state. The difference is computed through the norm 
%   ||A||=sum_kl |A_kl|^2. Obviously, the difference is zero
%   for Werner states. The form twirl(rho,d) makes it
%   possible to twirl a register of qudits with dimension d. 
%   Using the form twirl(rho,d,Nit) we can determine how many random
%   unitaries are used for twirling. (The algorithm
%   is not the straightforward integration of
%   the integral in the formula for twirling.)
%   The default value for Nit is 100.
%   For the algorithm see 
%   http://www.arxiv.org/abs/quant-ph/0609052.

function [r,difference]=twirl(rho,varargin);

if length(varargin)==0,
    % Dimension of quidits
    d=2;
    % Number of random unitaries used
    Nit=100;
elseif length(varargin)==1,
    d=varargin{1};
    Nit=100;
elseif length(varargin)==2,
    d=varargin{1};   
    Nit=varargin{2};
else
    error('Wrong number of input arguments');
end %if

x=[0 1;1 0];
z=[1 0;0 -1];
y=i*x*z;

[sy,sx]=size(rho);
N=log2(sx)/log2(d); 

U=zeros(d,d);
r=rho;
for n=1:Nit 
    % Initializing the random number generator based on the clock
    % is not a good idea, if the clock counts only seconds ...
    % Then after each restart we get back the same sequence of numbers
    % if the program does not run for more than a second.
    % Thus I omitted the following line:
    % if mod(n,100)==0,  rand('state',sum(100*clock));  end %if
    
    % Create a random dxd unitary
    % from d orthogonal vectors
    for k=1:d
        vv=randn(d,1)+i*randn(d,1); 
        for m=1:k-1
            vv=vv-U(:,m)*(U(:,m)'*vv);
        end %for
        U(:,k)=vv/sqrt(vv'*vv);
    end %for 
    
    UU=U;
    for n=2:N
        UU=kron(UU,U);   
    end %for
    r=(r+UU*r*UU')/2;
end %for

difference=trace((r-rho)*(r-rho)');


