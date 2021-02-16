% Test routines related to symmetric states, 
% especially sym2bipartite, sym2prodbasis and bipartite2prodbasis.
% Calculates the same thing in two ways and 
% stops if it finds an error

clear all
close all

% State (|01>+|10>)/sqrt(2), 1 qubit vs. 1 qubit

N=2;
M=1;
rho_sym=ketbra([0 1 0]); 
rho_sym
rho_bipartite=sym2bipartite(rho_sym,M);
rho_bipartite

rho_prodbasis=sym2prodbasis(rho_sym);
rho_prodbasis2=bipartite2prodbasis(rho_bipartite,M,N);

% Should be zero
difference=norm(rho_prodbasis-rho_prodbasis2)
if difference > 1e-15
    error(' ');
end

disp('------------------------')

% State (|1100>+permutations)/sqrt(6), 2 qubit vs. 2 qubits

N=4;
M=2;
rho_sym=ketbra([0 0 1 0 0]); 
rho_sym
rho_bipartite=sym2bipartite(rho_sym,M);
rho_bipartite

rho_prodbasis=sym2prodbasis(rho_sym);
rho_prodbasis2=bipartite2prodbasis(rho_bipartite,M,N);

% Should be zero
difference=norm(rho_prodbasis-rho_prodbasis2)
if difference > 1e-15
    error(' ');
end

disp('------------------------')

% State (|100>+|010>+|100>)/sqrt(3), 1 qubit vs. 2 qubits

N=3;
M=1;
rho_sym=ketbra([0 1 0 0]); 
rho_sym
rho_bipartite=sym2bipartite(rho_sym,M);
rho_bipartite

rho_prodbasis=sym2prodbasis(rho_sym);
rho_prodbasis2=bipartite2prodbasis(rho_bipartite,M,N);

% Should be zero
difference=norm(rho_prodbasis-rho_prodbasis2)
if difference > 1e-15
    error(' ');
end

disp('------------------------')

% State (|11>+|11>)/sqrt(2), 1 qubit vs. 1 qubit

N=2;
M=1;
rho_sym=ketbra([0 0 1]); 
rho_sym
rho_bipartite=sym2bipartite(rho_sym,M);
rho_bipartite

rho_prodbasis=sym2prodbasis(rho_sym);
rho_prodbasis2=bipartite2prodbasis(rho_bipartite,M,N);

% Should be zero
difference=norm(rho_prodbasis-rho_prodbasis2)
if difference > 1e-15
    error(' ');
end

disp('------------------------')

% State (|100>+|010>+|100>)/sqrt(3), 2 qubits vs. 1 qubits

N=3;
M=2;
rho_sym=ketbra([0 1 0 0]); 
rho_sym
rho_bipartite=sym2bipartite(rho_sym,M);
rho_bipartite

rho_prodbasis=sym2prodbasis(rho_sym);
rho_prodbasis2=bipartite2prodbasis(rho_bipartite,M,N);

% Should be zero
difference=norm(rho_prodbasis-rho_prodbasis2)
if difference > 1e-15
    error(' ');
end

disp('------------------------')

% Random N-qubit state

N=5;
M=3;
rho_sym=rdmat(1,N+1); 
rho_sym
rho_bipartite=sym2bipartite(rho_sym,M);
rho_bipartite

rho_prodbasis=sym2prodbasis(rho_sym);
rho_prodbasis2=bipartite2prodbasis(rho_bipartite,M,N);

% Should be zero
difference=norm(rho_prodbasis-rho_prodbasis2)
if difference > 1e-15
    error(' ');
end



