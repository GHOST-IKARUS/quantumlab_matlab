% BES_private   Bound entangled state based on private states
%    The state is for a DxD system and is given in
%    https://arxiv.org/abs/2002.12409
%    https://arxiv.org/abs/1309.7992
%    The function must be used as BES_private(D), where
%    the sate will be of dimension DxD. 
%
%    BES_private(n,v) with v=0 provides the state
%    corresponding to the order of subsystems
%    ABA'B', which is also used in the paper.
%    BES_private(n,v) with v=1 provides the state
%    corresponding to the order of subsystems
%    AA'BB', which is better to study entanglement between AA' and BB'.
%    If the parameter v is omitted then its value is 
%    taken to be 0.
%
%    See also example_BES_private.

function rho=BES_private(D,varargin)

if isempty(varargin),
    % Default order
    oder_of_subsystems=0;
else
    oder_of_subsystems=varargin{1};
end %if

d=D/2;

% U is a dxd matrix
% Quantum Fourier transform
i=sqrt(-1);
for j=1:d
    for k=1:d
        U(k,j)=sqrt(1/d)*exp(i*2*pi*j*k/d);
    end
end

p1=sqrt(d)/(1+sqrt(d));
p2=1-p1;

sum1=0;
for i=0:d-1
    for j=0:d-1
        v1=mkron(uvec(1,2),uvec(1,2),uvec(i+1,d),uvec(j+1,d));
        v2=mkron(uvec(2,2),uvec(2,2),uvec(i+1,d),uvec(j+1,d));
        sum1=sum1+ketbra(v1)+ketbra(v2);
    end
end

sum_u=0;
for i=0:d-1
    for j=0:d-1
        v1=mkron(uvec(1,2),uvec(1,2),uvec(i+1,d),uvec(j+1,d));
        v1b=mkron(uvec(2,2),uvec(2,2),uvec(j+1,d),uvec(i+1,d));
        v2=mkron(uvec(2,2),uvec(2,2),uvec(j+1,d),uvec(i+1,d));
        v2b=mkron(uvec(1,2),uvec(1,2),uvec(i+1,d),uvec(j+1,d));
        sum_u=sum_u+U(i+1,j+1)*(ket(v1)*bra(v1b))+...
              conj(U(i+1,j+1))*(ket(v2)*bra(v2b));
    end
end

sum2=0;
for i=0:d-1
    v1=mkron(uvec(1,2),uvec(2,2),uvec(i+1,d),uvec(i+1,d));
    v2=mkron(uvec(2,2),uvec(1,2),uvec(i+1,d),uvec(i+1,d));
    sum2=sum2+ketbra(v1)+ketbra(v2);
end

sum3=0;
for i=0:d-1
    for j=0:d-1
        v1=mkron(uvec(1,2),uvec(2,2),uvec(i+1,d),uvec(i+1,d));
        v1b=mkron(uvec(2,2),uvec(1,2),uvec(j+1,d),uvec(j+1,d));
        v2=mkron(uvec(2,2),uvec(1,2),uvec(j+1,d),uvec(j+1,d));
        v2b=mkron(uvec(1,2),uvec(2,2),uvec(i+1,d),uvec(i+1,d));
        sum3=sum3+U(i+1,j+1)*(ket(v1)*bra(v1b))+...
             conj(U(i+1,j+1))*(ket(v2)*bra(v2b));
    end
end


rho=p1/2/d^2*sum1+p1/2/d/sqrt(d)*sum_u+p2/2/d*sum2+p2/2/d*sum3;

if oder_of_subsystems==1

   rho=reorder(rho,[4 2 3 1],[2 2 d d]);

end
