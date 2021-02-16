% elin    Convex roof of the linear entropy (linear entropy of
%         entanglement)
%         See https://arxiv.org/abs/1409.3806
%         For an alternative realization see
%         http://de.mathworks.com/matlabcentral/fileexchange/47823-corona-convex-roof-numerical-analysis?requestedDomain=www.mathworks.com

function E_lin=elin(rho_full)

% TEST 1=ON,0=OFF
TEST=1;

% Treshold to consider an eigenvalue zero
lambdatreshold=0.000001;

[sy,sx]=size(rho_full);

% Dimension of qudits
dq=floor(sqrt(sy)+0.5);

% Number of two-qudit units
% for the PPT symmetric extension
% N=2, first nontrivial level
N=2;

% Calculate the matrix for dimension reduction
% and the density matrix in the new basis
if 1 %rank(rho_full)<dq^2
    [v,dd]=eig(rho_full);
    dd=diag(dd);
    % Tereshold to look for zero eigenvalue
    index=find(abs(dd)>lambdatreshold);
    
    Ureduce=[];
    for n=1:length(index)
        Ureduce=[Ureduce,v(:,index(n))];
    end %for
    
    % d=rank of the state
    d=length(index);
    rho=Ureduce'*rho_full*Ureduce;
    
    % Just in case, normalize rho
    rho=diag(diag(rho));
    rho=real(rho);
    rho=nm(rho);
       
    %Ur=runitary(1,d);
    %Ureduce=Ureduce*Ur;
    
    %rho=diag(dd(index));
    % Note: rho_full=Ureduce*rho*Ureduce'
else
    d=dq^2;
    rho=rho_full;
    Ureduce=eye(d);
end %if

% Just to be sure, normalize it
rho=nm(rho);

% CHOOSE SOLVER
%yalmip('solver','sdpt3')
%yalmip('solver','sedumi')
sdpsettings('verbose',1);

% Observable in the expectation value
Obs=zeros(dq^4,dq^4);
oo=sud(dq);
for n=1:dq^2-1
    Obs=Obs+mkron(oo(:,:,n),eye(dq),oo(:,:,n),eye(dq));
end %for

% ONLY FOR TEST
% Compute sum_k G_k^2
%sGk2=0;
%for n=1:dq^2-1
%    sGk2=sGk2+oo(:,:,n)^2;
%end %for
%trace_sGk2=trace(sGk2)

% Reduce to the range of rho
Obs0=Obs;
Obs=kron(Ureduce,Ureduce)'*Obs*kron(Ureduce,Ureduce);
Obs=(Obs+Obs')/2;

if d==2
    rhols=sdpvar(N+1,N+1,'hermitian','complex');
    %rhol=0*sdpvar(d^N,d^N,'hermitian','complex');
    rhol=zeros(d^N,d^N);
    for n=0:N
        for m=0:N
            rhol=rhol+rhols(n+1,m+1)*dstate(n,N)*dstate(m,N)';
        end %for
    end %for
elseif (d==3 || d==4 || d==9 || d==6) || N==2
    
    Ps=proj_sym(N,d);
    rPs=rank(Ps);
    rhols=sdpvar(rPs,rPs,'hermitian','complex');
    [v,dd]=eig(Ps);
    index=find(diag(dd)>0.99);
    v=v(:,index);
        
    for n=1:rPs
        for m=1:rPs
            if TEST
                disp(['General state-symmetric state ' num2str(n) ' ' num2str(m)]);
            end %if
            if n==1 && m==1
                rhol=rhols(1,1)*v(:,1)*v(:,1)'; %Intialization
            else
                rhol=rhol+rhols(n,m)*v(:,n)*v(:,m)';
            end %for
            
        end %for
    end %for
    
else
    
    error(['Error with d and N: the procedure is not yet implemented for these paramaters: d=' num2str(d) '.'])
    
end %if

if TEST
    disp(['1: general state-symmetric state conversion is done'])
end %if

% Computing the reduced states
[sx,sy]=size(rhol);
N=log2(sx)/log2(d);
N=floor(N+0.5);
listneg=1:2;
Nr=length(listneg);
if N==2
    rho_twoqudit=rhol;
else
    rho_twoqudit=sdpvar(d,d,'hermitian','complex');
    % TEST
    %[sy,sx]=size(rhol);
    %R=randn(sx,sy);
    %rho2=R'*rhol*R;
    %disp(['11'])
    for k=0:d^Nr-1
        for l=0:d^Nr-1
            if TEST
                disp(['Two-body ' num2str(k) ' ' num2str(l)]);
            end %if
            rr=0*rhol(1,1);
            i1=1+k*d^(N-Nr);
            i2=1+l*d^(N-Nr);
            for n=0:d^(N-Nr)-1
                rr=rr+rhol(i1+n,i2+n);
            end %for
            rho_twoqudit(k+1,l+1)=rr;
        end %for
    end %for
end %if

if TEST
    disp(['2: Computing two-body reduced state is done.'])
end %if

if N==2
    rho_singlequdit=[];
    for n=1:d
        rho_singlequdit_row=[];
        for m=1:d
            rho_singlequdit_row=[rho_singlequdit_row,trace(rhol((n-1)*d+1:(n-1)*d+d,(m-1)*d+1:(m-1)*d+d))];
        end %for
        rho_singlequdit=[rho_singlequdit;rho_singlequdit_row];
    end %for
else
    listneg=1;
    Nr=length(listneg);
    rho_singlequdit=sdpvar(d,d,'hermitian','complex')
    for k=0:d^Nr-1
        for l=0:d^Nr-1
            if TEST
                disp(['Single-particle reduced state ' num2str(k) ' ' num2str(l)]);
            end %if
            rr=0*rhol(1,1);
            i1=1+k*d^(N-Nr);
            i2=1+l*d^(N-Nr);
            for n=0:d^(N-Nr)-1
                rr=rr+rhol(i1+n,i2+n);
            end %for
            rho_singlequdit(k+1,l+1)=rr;
        end %for
    end %for
end %if

if TEST
    disp(['3: Computing single-particle reduced state is done.'])
end %if

if N==2
    rhol_PT=[];
    for n=1:d
        rhol_PT_row=[];
        for m=1:d
            rhol_PT_row=[rhol_PT_row,rhol((n-1)*d+1:(n-1)*d+d,(m-1)*d+1:(m-1)*d+d).'];
        end %for
        rhol_PT=[rhol_PT;rhol_PT_row];
    end %for
else
    % Partial transpose
    rhol_PT=rhol; % NOT VERY NICE, BUT WORKS.
    if N==2
        N1=1;
        N2=1;
    elseif N==3
        N1=2;
        N2=1;
    elseif N==4
        N1=1;
        N2=3;
    end %if
    AM=kron(1:d^N,ones(d^N,1));
    BM=kron(ones(1,d^N),(1:d^N).');
    AT=pt_nonorm(AM,1:N1,d);
    BT=pt_nonorm(BM,1:N1,d);
    for m=1:d^N
        for n=1:d^N
            if TEST
                disp(['Partial transpose ' num2str(m) ' ' num2str(n)]);
            end %if
            rhol_PT(m,n)=rhol(AT(m,n),BT(m,n));
        end %for
    end %for
end %if

if TEST
    disp(['4: Comptuing the partial transpose is done'])
end %if

% Of course, we cannot juts write F=F+set(rho_singlequdit==rho);
F=set(rhols>=0)+set(trace(rhols)==1);
for n=1:d
    for m=1:n
        if TEST
           disp(['Constraints ' num2str(n) ' ' num2str(m)]);
       end %if
        F=F+set(rho_singlequdit(m,n)==rho(m,n));
    end %for
end %for

F=F+set(rhol_PT>=0);

if TEST
    disp(['5: Contstaints are done'])
end %if

diagnostic=solvesdp(F,-real(trace(Obs*rho_twoqudit)));

% if numerics had problems, give back NaN
if diagnostic.problem~=0
    E_Lin=NaN;
end %if

supremum=double(trace(Obs*rho_twoqudit));

E_lin=1-1/dq-(1/2)*supremum;

% Only for test
if TEST
    
    diagnostic=diagnostic
    
    supremum=supremum
    
    supremum2_should_be_equal_to_supremum=double(trace(Obs0*kron(Ureduce,Ureduce)*rho_twoqudit*kron(Ureduce,Ureduce)'))
    
    E_lin=E_lin
    
    r2=double(rho_twoqudit);
    
    r2=kron(Ureduce,Ureduce)*r2*kron(Ureduce,Ureduce)';
    
    % A single copy of the bipartite state
    r1=keep(r2,2:-1:1,dq);
    
    difference_should_be_zero=max(max(abs(r1-rho_full)))
    
    neg_rho=negativity(rho_full,1,dq)
    
end %if

if E_lin<-1e-5
    error('Error: Negative enanglement!')
end %if



