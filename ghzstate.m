%ghzstate   Defines a GHZ state.
%   ghzstate(n,d) gives the state vector for an n-qubit GHZ state of
%   d dimensional particles. If d is omitted, then it is taked to be 2.
%   If argument n is omitted than the default is taken to be
%   the value of global variable N.

function g=ghzstate(varargin)

if isempty(varargin)
    global N;
    
elseif length(varargin)==1
    % GHZ state of qubits
    N=varargin{1};
    g=zeros(2^N,1);
    g(1)=1/sqrt(2);
    g(end)=1/sqrt(2);
    
elseif length(varargin)==2
    % Higher dimensional GHZ state
    N=varargin{1};
    d=varargin{2};
    g=zeros(d^N,1);
    % Compute 1+d+d^2+...+d^(N-1)
    step=sum(d.^(0:N-1));
    g(1+(0:d-1)*step)=1/sqrt(d);
else
    error('Wrong number of input arguments.')
end %if

