% interact   Interaction between two qudits
%   interact (op1,op2,n1,n2,n) gives an operator
%   acting on qudits n1 and n2, respectively, with
%   operators op1 and op2. n is the number of qudits.
%   If argument n is omitted than the default
%   is taken to be the value of global variable N.
%   The dimension of the qudit is obtained from the size of op1.
%   (For qubits op1 is 2x2, for qutrits it is 3x3, etc.)

function op=interact(op1,op2,n1,n2,varargin)

if isempty(varargin),
    global N;
else
    if length(varargin)~=1,
        error('Wrong number of input arguments.')
    end %if
    N=varargin{1};
end %if

d=max(size(op1));

m1=min(n1,n2);
m2=max(n1,n2);
op=mkron(eye(d^(N-m2)),op2,eye(d^(m2-m1-1)),op1,eye(d^(m1-1)));
