% fisherwit_spinsq   Plots the optimal lower bounds on the quantum Fisher
%                   information F_Q[rho,Jy] as the function of
%                   the expectation values of Jz and Jx^2. This way we can witness
%                   the quantum Fisher information based
%                   on measuring some expectation values.
%                   The programs plot the curve for a set of
%                   particle numbers.
%
%                   The program was used for the publication
%
%                   Apellaniz, Kleinmann, Guhne, Toth,
%                   Optimal detection of metrologically useful entanglement
%                   http://arxiv.org/abs/1511.05203

tic

format compact


% General states vs. symmetric states (faster with symmetric states)
GeneralStates=0;

% Sparse matrices ON/OFF (faster with ON)
SPARSEMATRICES_ON=0;

% Parameters for spin squeezing
N_range=[40 60 80 100]; % range of particle numbers
rstep_initial=[1 1]; % Initial step size
rstep_divide=2; % Step size is divided by this value, when it is decreased
numberofunsuccesfultrials=[80 10 10 10 10 10 10 10 10 10 10 10 10]; % After this many unsuccesful trials, it decreases the stepsize
rstep_final=0.01; % The step size is decreased until this limit
mu_resolution=200;50;% Number of points considered between the two extreme value for mu
Warray_str='{Jx^2,Jz}'; % Witnesses
scalingarray_str='[N,N]';

% Spin squeezing parameter
xi2=0.1514;

% Basic relations
%
% xi2=N var(Jx)/<Jz>^2.
% var(Jx)=xi2/N<Jz>^2
% <Jz>=alpha N/2
% var(Jx)=xi2*N*alpha^2/4

alpha=0.85;

warray0=[xi2/4*alpha^2,1/2*alpha]; % Expectation values, normalized

% Check physicality
NN=min(N_range);
j=NN/2;
minvarJz=Fj(warray0(2)*NN/j,j)*j
varJz=warray0(1)*NN
if varJz<minvarJz
    error('Aphysical spin squeezed state')
end %if

% Display the Pezze-Smerzi bound
FQ_Pezze_Smerzi=(warray0(2))^2/warray0(1)

FQmin_optim_array=[];
rvector_array=[];
rvector=zeros(1,length(warray0));

for iN=1:length(N_range)
    
    N=N_range(iN)
    
    mu_range=-N/2:N/(mu_resolution+1):N/2; % range for mu
    
    if GeneralStates==1
        % General states
        paulixyz;
        Jx=coll(x,N)/2;
        Jy=coll(y,N)/2;
        Jz=coll(z,N)/2;
        A=coll(y,N)/2;
        Identity=coll(e,N);
    else
        % Symmetric states
        d=N+1;
        [Jx,Jy,Jz]=su2(d);
        A=Jy;
        Identity=eye(d,d);
    end %if
    
    if SPARSEMATRICES_ON
        Jx=sparse(Jx);
        Jy=sparse(Jy);
        Jz=sparse(Jz);
        Identity=sparse(Identity);
    end %if
    
    % INPUT
    % Witnesses and values
    Warray=eval(Warray_str);
    warray=warray0.*eval(scalingarray_str);
    
    rstep=rstep_initial;
    
    FQmin_optim=-Inf;
    
    trialcounter=0;
    
    % Looking for the optimal r,
    % exploiting the concavity, convexity
    % properties.
    while max(rstep)>rstep_final
        
        r0=rvector;
        Wrsum=0;
        wrsum=0;
        for ri2=1:length(Warray)
            Wrsum=Wrsum+r0(ri2)*Warray{ri2};
            wrsum=wrsum+r0(ri2)*warray(ri2);
        end %for
        
        % Compute the derivative in r at r0
        % Value for r0
        flag_degeneracy=1;
        r00=r0;
        while flag_degeneracy
            r=r0;
            %Wr=r*W;
            witA=[];
            op_muindependent=Wrsum-4*A^2;
            % Loop for different mu values
            flag_degeneracy=0;
            for Mu=mu_range
                H=op_muindependent+8*Mu*A;
                if SPARSEMATRICES_ON
                    witA=[witA,max(real(eigs(H,1,'LR')))-4*Mu^2];
                else
                    witA=[witA,max(real(eig(H)))-4*Mu^2];
                end %if
            end %for
        end %while
        
        % Store the witA and wrsum corresponding to the optimum
        witA0=witA;
        wrsum0=wrsum;
        
        maxwitA=max(witA);
        FQmin0=wrsum-maxwitA;
        
        % Value for r0+deltar
        deltar=rstep.*randn(1,length(Warray));
        
        r=r0+deltar;
        Wrsum=0;
        wrsum=0;
        for ri2=1:length(Warray)
            Wrsum=Wrsum+r(ri2)*Warray{ri2};
            wrsum=wrsum+r(ri2)*warray(ri2);
        end %for
        
        witA=[];
        op_muindependent=Wrsum-4*A^2;
        % Loop for different mu values
        flag_degeneracy=0;
        for Mu=mu_range
            H=op_muindependent+8*Mu*A;
            if SPARSEMATRICES_ON
                witA=[witA,max(real(eigs(H,1,'LR')))-4*Mu^2];
            else
                witA=[witA,max(real(eig(H)))-4*Mu^2];
            end %if
        end %for
               
        maxwitA=max(witA);
        FQmin1=wrsum-maxwitA;
        
        % Display
        disp(['r0=' num2str(r0) ' / r=' num2str(r) ' / FQ0=' num2str(FQmin0) ' FQ=' num2str(FQmin1) ' / deltar=' num2str(rstep) ]);
        
        if FQmin1>FQmin0
            rvector=r;
            trialcounter=0;
        else
            trialcounter=trialcounter+1;
            if trialcounter>numberofunsuccesfultrials(iN),
                rstep=rstep/rstep_divide;
                trialcounter=0;
            end %if
        end %if
                
    end %while
    
    rvector
    
    FQmin_optim_array=[FQmin_optim_array,FQmin0];
    
    rvector_array=[rvector_array,rvector.'];
    
end %for

p=plot(N_range,FQmin_optim_array./N_range,'o-')

t1=xlabel('$N$')
t2=ylabel('$F_Q[\varrho,J_y]/N$')
%axis([0 1 0 1])
set(p,'LineWidth',2);
set(gca,'LineWidth',2);
set(gca,'FontSize',20);
set(t1,'Interpreter','Latex');
set(t2,'Interpreter','Latex');
set(t1,'FontSize',24);
set(t2,'FontSize',24);

toc


