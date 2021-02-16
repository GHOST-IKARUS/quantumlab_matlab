% BES_metro   Bound entangled state for robust metrology
%    The state is for a DxD system and is given in
%    https://arxiv.org/abs/2002.12409
%    The function must be used as BES_metro2020(n), where
%    the sate will be of dimension 2^(n+1)x2^(n+1).
%    Works for D=6 or power of two as 4,8,16,32,64,128,etc.
%
%    BES_metro2020(n,v) with v=0 provides the state
%    corresponding to the order of subsystems
%    ABA'B', which is also used in the paper.
%    BES_metro2020(n,v) with v=1 provides the state
%    corresponding to the order of subsystems
%    AA'BB', which is better to study entanglement between AA' and BB'.
%    If the parameter v is omitted then its value is
%    taken to be 0.
%
%    See also example_BES_metro.

function rho=BES_metro(D,varargin)

if isempty(varargin),
    % Default order
    oder_of_subsystems=0;
else
    oder_of_subsystems=varargin{1};
end %if

if D==6
    
    d=3;
    
    phi0=0;
    Q=zeros(3,3,3);
    
    k=1;
    f=2*pi/3*k+phi0;
    Q(1,:,:)=[cos(f) sin(f) 0; sin(f) -cos(f) 0;0 0 1];
    
    k=2;
    f=2*pi/3*k+phi0;
    Q(2,:,:)=[cos(f) sin(f) 0; sin(f) -cos(f) 0;0 0 1];
    
    k=3;
    f=2*pi/3*k+phi0;
    Q(3,:,:)=[cos(f) sin(f) 0; sin(f) -cos(f) 0;0 0 1];
    
    p1=sqrt(d)/(1+sqrt(d));
    p2=1-p1;
    
    sum1=0;
    for i=0:d-1
        for j=0:d-1
            zij=mkron(uvec(1,2),uvec(1,2),uvec(i+1,d),uvec(j+1,d))/sqrt(2);
            for k=0:d-1
                % ik or ki??
                zij=zij+Q(j+1,i+1,k+1)*mkron(uvec(2,2),uvec(2,2),uvec(j+1,d),uvec(k+1,d))/sqrt(2);
            end
            sum1=sum1+ketbra(zij);
        end
    end
    
    s0=(mkron(uvec(1,2),uvec(2,2),uvec(1,d),uvec(1,d))...
        +mkron(uvec(1,2),uvec(2,2),uvec(2,d),uvec(2,d)))/sqrt(2);
    s1=(mkron(uvec(1,2),uvec(2,2),uvec(1,d),uvec(2,d))...
        -mkron(uvec(1,2),uvec(2,2),uvec(2,d),uvec(1,d)))/sqrt(2);
    s2= mkron(uvec(1,2),uvec(2,2),uvec(3,d),uvec(3,d))/sqrt(2);
    sum2=ketbra(s0)+ketbra(s1)+ketbra(s2);
    
    sum3=0;
    for i=0:d-1
        sum3=sum3+ketbra(mkron(uvec(2,2),uvec(1,2),uvec(i+1,d),uvec(i+1,d)));
    end
    
    rho=p1/d^2*sum1+p2/2/d*sum2+p2/2/d*sum3;
    
    if oder_of_subsystems==1
        
        rho=reorder(rho,[4 2 3 1],[2 2 d d]);
        
    end
    
    
elseif D>=4 && round(log2(D))==log2(D)
    
    n=round(log2(D/2));
    d=D/2;
    
    
    P=zeros(2^n,2^n,2^n);
    X=[0 1;1 0];
    for k=0:2^n-1
        b=dec2bin(k,n);
        p=str2num(b(1));
        Pk=X^p;
        for m=2:n
            p=str2num(b(m));
            Pk=kron(Pk,X^p);
        end
        P(:,:,k+1)=Pk;
    end
    Q=P;
    
    % Obtain t_k for k=0:d-1
    tk=zeros(d^2,d);
    E=eye(2^n,2^n);
    for k=0:d-1
        t=zeros(d^2,1);
        for i=1:2^n
            for j=1:2^n
                t=t+P(i,j,k+1)*kron(uvec(i,d),uvec(j,d))/sqrt(d);
            end
        end
        tk(:,k+1)=t;
    end
    sk=tk;
    
    p1=sqrt(d)/(1+sqrt(d));
    p2=1-p1;
    
    sum1=0;
    for i=0:d-1
        for j=0:d-1
            zij=mkron(uvec(1,2),uvec(1,2),uvec(i+1,d),uvec(j+1,d))/sqrt(2);
            for k=0:d-1
                % ik or ki??
                zij=zij+Q(j+1,i+1,k+1)*mkron(uvec(2,2),uvec(2,2),uvec(j+1,d),uvec(k+1,d))/sqrt(2);
            end
            sum1=sum1+ketbra(zij);
        end
    end
    
    % Is it 0 1 0 0 or 0 0 1 0?
    sum2=0;
    for k=0:d-1
        sum2=sum2+kron(ketbra([0 1 0 0]),ketbra(sk(:,k+1)));
    end
    
    sum3=0;
    for i=0:d-1
        sum3=sum3+ketbra(mkron(uvec(2,2),uvec(1,2),uvec(i+1,d),uvec(i+1,d)));
    end
    
    rho=p1/d^2*sum1+p2/2/d*sum2+p2/2/d*sum3;
    
    if oder_of_subsystems==1
        
        rho=reorder(rho,[4 2 3 1],[2 2 d d]);
        
    end
    
else
    
    error('D must be 6 or power of two as 4,8,16,32,64,128,etc.')
    
end




