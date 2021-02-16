% fisherwit_dicke   Plots the lower bounding on the quantum Fisher 
%                   information F_Q[rho,Jy] as the function of 
%                   the fidelity with respect to symmetric 
%                   Dicke states. This way we can witness
%                   the quantum Fisher information based 
%                   on measuring the fidelity.
%
%                   The program plots the curve for N=6 and 12
%                   particles.
%
%                   The program was used for the publication 
%
%                   Apellaniz, Kleinmann, Guhne, Toth, 
%                   Optimal detection of metrologically useful entanglement
%                   http://arxiv.org/abs/1511.05203

% Geza Toth (2015)

clear all
close all

tic

format compact

% Parameters for fast calculattion
% If you increase them, you get more exact figures,
% but the program will be slower.
% For larger N, r_range_max must be larger.
N_range=[6 12];  % range of particle numbers
r_range_min=0; % r plays the role of the slope (only minimum, maximum counts)
r_range_max=1000; % upper bound for r
deltar=0.01; % For computing the derivative
r_resolution=0.1; % For the resolution reuqired
mu_resolution=300; % Number of points considered between the two extreme value for mu

% Range of fidelities for which we carry out the calculation
fidelity_range=[0:0.02:0.8,0.81:0.01:0.9,0.91:0.005:1];

for N=N_range
    
    mu_range=-N/2:N/2/(mu_resolution+1):N/2; % range for mu
    
    FQmin_optim_array=[];
    
    GeneralStates=0;
    if GeneralStates==1
        % General states
        paulixyz;
        A=coll(x,N)/2;
        Identity=coll(e,N)/2;
        Pdicke=ketbra(dstate(N/2,N));
    else
        % Symmetric states
        d=N+1;
        [Jx,Jy,Jz]=su2(d);
        A=Jx;
        Identity=eye(d,d);
        Pdicke=zeros(N+1,N+1);
        Pdicke(N/2+1,N/2+1)=1; %Dicke state
    end %if
    
    N
    
    % Larger Fidelity allowing still FQ=0.
    fy=grstate(-Jy);
    MaxFidelity_Allowing_FQ0=trace(ketbra(fy)*Pdicke)
    
    % Loop for different exp. values of the projector
    for fidelity=fidelity_range
        
        %FQmin_optim=-Inf;
        
        rmin=r_range_min;
        rmax=r_range_max;
        
        % Looking for the optimcal r,
        % exploiting the concavity, convexity
        % properties.
        % Interval halving method
        while abs(rmin-rmax)>r_resolution
            
            r0=(rmin+rmax)/2;
            
            % Compute the derivative in r at r0
            % Value for r0
            r=r0;
            Wr=r*Pdicke;
            witA=[];
            % Loop for different mu values
            for Mu=mu_range
                H=Wr-4*A^2+8*Mu*A-4*Mu^2*Identity;
                [v,d]=eig(H);
                [values,index]=max(real(diag(d)));
                witA=[witA,values(1)];
                % Degeneracy check is not needed if
                % GeneralStates==0
                if length(index)>1
                    error('Degeneracy!')
                end %if
            end %for
            maxwitA=max(witA);
            FQmin0=r*fidelity-maxwitA;
            % Value for r0+deltar
            r=r0+deltar;
            Wr=r*Pdicke;
            witA=[];
            % Loop for different mu values
            for Mu=mu_range
                H=Wr-4*A^2+8*Mu*A-4*Mu^2*Identity;
                [v,d]=eig(H);
                [values,index]=max(real(diag(d)));
                witA=[witA,values(1)];
                % Degeneracy check is not needed if
                % GeneralStates==0
                if length(index)>1
                    error('Degeneracy!')
                end %if
            end %for
            maxwitA=max(witA);
            FQmin1=r*fidelity-maxwitA;
            
            deriv=(FQmin1-FQmin0)/deltar;
            
            if deriv>0
                rmin=r;
            else
                rmax=r;
            end %if
            
        end %while
        
        FQmin_optim_array=[FQmin_optim_array,FQmin0];
        
    end %for
    
    h=plot(fidelity_range,FQmin_optim_array/N^2);
    set(h,'LineWidth',2);
    
    hold on
    
end %for

grid

fs=30; gfs=26;
t1=xlabel('$F_{\rm Dicke}$');
t2=ylabel('$F_Q[\varrho,J_y]/N^2$');
axis([0 1 0 1])
set(gca,'LineWidth',2);
set(gca,'FontSize',gfs);
set(t1,'Interpreter','Latex');
set(t2,'Interpreter','Latex');
set(t1,'FontSize',fs);
set(t2,'FontSize',fs);

% Move figure up to make to make space for label
p=get(gca,'Position');
set(gca,'Position',p+[0 0.05 0 0]);

toc

