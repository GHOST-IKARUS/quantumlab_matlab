% MATLAB routines for quantum mechanics.
% 
% General
%  contents              - List of commands
%  ver                   - Version
%  example1              - Example with two-state systems
%  example2              - Example with multi-qubit states
%  example3              - Example about operators in spin chains
%
%
% Bras and kets
%  ket                   - Creating a normalized column vector (Dirac's ket)
%  bra                   - Creating a normalized row vector    (Dirac's bra)
%  ketbra                - Creating a projector from a vector
%  ketbra2               - Creating a projector from the argument 
%                          if the argument is a vector.
%  braket                - Dirac's braket
%  ex                    - Expectation value
%  va                    - Variance
%  nm                    - Normalization
%  uvec                  - Vector with a single 1 element
%
%
% Reordering a qudit register
%  mkron                 - Kronecker product with several arguments
%  pkron                 - Multiplying a matrix with itself given times
%                          using the Kronecker product
%  remove                - Reduced density matrix in a qudit register
%  keep                  - Reduced density matrix in a qudit register
%  reorder               - Reorder a register of qudits
%  reordermat            - The operator corresponding to reordering
%                          a register of qudits
%  reordervec            - Transformation vector for reordering 
%                          a register of qudits
%  shiftquditsleft       - Shift qudits to the left
%  shiftquditsright      - Shift qudits to the right
%  swapqudits            - Swap two qudits of a quantum state
%
%
% Interesting quantum states, quantum gates and operators
%  ghzstate              - Greenberger-Horne-Zeilinger state
%  cstate                - Cluster state
%  rstate                - Ring cluster state
%  gstate                - Graph state
%  gstate_stabilizer     - Stabilizer of a graph state
%  wstate                - W state
%  dstate                - Symmetric Dicke state
%  singlet               - Singlet state
%  mestate               - Maximally entangled state
%  mmstate               - Density matrix for the maximally mixed state
%  smolinstate           - Smolin's four-qubit bound entangled state
%  BES_Horodecki3x3      - Horodecki's 3x3 bound entangled state
%  BES_Horodecki2x4      - Horodecki's 2x4 bound entangled state
%  BES_private           - Bound entangled states based on private states.
%                          See example_BES_private.
%  BES_metro4x4          - 4x4 Bound entangled state for quantum metrology
%                          See example_BES_metro4x4.
%  BES_metro             - Family of bound entangled states for quantum metrology
%                          See example_BES_metro.
%  BES_UPB3x3            - 3x3 UPB bound entangled state
%  BES_Breuer            - Breuer's bound entangled state
%  BES_Watrous           - Watrous' bound entangled state
%  U_CNOT                - 4x4 unitary matrix of a CNOT gate
%  U_H                   - 2x2 unitary matrix for the Hadamard gate
%  paulixyz              - Define Pauli matrices x,y,z and e=eye(2)
%  Jxyz                  - Collective angular momentum operators for N
%                          qubits
%  su2                   - SU(2) generators for matrices of a given size
%  su3                   - Define the SU(3) generators m1,...,m8 and ee=eye(3)
%  su3_alternative       - Define alternative SU(3) generators
%  sud                   - Define the SU(d) generators
%
%
% Continuous variable (CV) systems 
%  aop                   - annihilation operator
%  nop                   - number operator
%  xop                   - x operator
%  pop                   - p operator
%  cohstate              - coherent state
%
%
% Formatted input/output
%  printv                - Print state vector in product basis
%  decompose             - Display pauli decomposition of a matrix
%  paulistr              - Convert symbolic string to operator 
%
%
% Two-qudit interactions and spin chains
%  quditop               - Operator acting on a qudit of a qudit register
%  twoquditop            - Operator acting on two qudits of a qudit register
%  coll                  - Define a collective multi-qudit operator
%  interact              - Two-qudit interaction acting on given qudits
%  nnchain               - Spin chain Hamiltonian
%  nnchainp              - Spin chain Hamiltonian with a periodic boundary
%                          condition
%  ising                 - Ising spin chain Hamiltonian
%  isingp                - Ising spin chain Hamiltonian with a periodic 
%                          boundary condition
%  ising_ground          - Ground state energy of Ising model
%  ising_free            - Free energy in thermal state
%  ising_thermal         - Internal energy in thermal state
%  ising_classical_ground - Ground state energy for the classical Ising model
%  heisenberg            - Heisenberg spin chain Hamiltonian
%  heisenbergp           - Heisenberg spin chain Hamiltonian with a periodic 
%                          boundary condition
%  xy_classical_ground   - Ground state energy for the classical xy model
%  orthogobs             - Orthogonal observables for a qudit
%
%
% Entanglement
%  pt                    - Partial transpose for a qudit register
%  pt_nonorm             - Like ppt, but without normalization. 
%                          Works with the sdp frontend yalmip.
%  negativity            - Compute the negativity of the density
%                          matrix
%  concurrence           - Concurrence for a two-qubit matrix
%  realign               - Realignment of a density matrix
%  ccnr                  - Computable Cross Norm - Realignment criterion
%  mrealign              - Realignment for multiqudit systems
%  schmidt               - Schmidt coefficients for a pure state
%  maxsep                - Maximum of an operator for separable states
%  maxsymsep             - The same as maxsep but for
%                          permutationally invariant sep. states
%  maxbisep              - Maximum of an operator for biseparable states
%  maxb                  - Like maxbisep, but for all bipartitions
%  overlapb              - Maximum overlap of a pure state with bisep. states
%  maxppt                - Maximum of an operator expectation value 
%                          for states with a positive partial transpose
%  example_maxppt        - Example for the usage of maxppt
%  optwitness            - Obtaining optimal entanglement witnesses
%  example_optwitness    - Example for the usage of optwitness
%  elin                  - Linear entropy of entanglement using
%                          semidefinite programing. Works faster for low 
%                          rank states.
%  example_elin          - Example for the elin function
%  example_elin2         - Second example for the elin function
%
%
% Spin squeezing
%  Fj                    - The function Fj(x) from Sorensen-Molmer
%                          extreme spin squeezing paper
%  Fj_inv                - Inverse of Fj
%  Fj_approx             - A lower bound on Fj(x), which is fast to compute
%  example_Fj            - Example for the Fj function
%  example_Fj_inv        - Example for the Fj_inv function
%  example_Fj_approx     - Example for the Fj_approx function
%  optspinsq             - Optimal spin squeezing inequalities
%
%
% Quantum metrology
%  fisher                - Quantum Fisher information (qFi) 
%  sld                   - Symmetric logarithmic derivative (SLD)
%  skewinf               - Wigner-Yanase skew information
%  fisherwit_dicke       - Bounding the qFi from below based on the fidelity
%                          with respect to Dicke states
%  fisherwit_spinsq      - Bounding the qFi from below based on spin
%                          squeezing experiments
%
%
% Random vectors, matrices and operations
%  rvec                  - Random state vector for a given number of
%                          qudits
%  rproduct              - Random product state vector for a given
%                          number of qudits
%  rdmat                 - Random density matrix for a qudit register
%  runitary              - Random unitary for a qudit register
%  rhermitian            - Random Hermitian matrix
%  twirl                 - Twirling
%  twirl2                - How close is a state to Werner states
%
%
% Miscellaneous simple commands
%  proj_sym              - Projector to the symmetric subspace
%  proj_asym             - Projector to the antisymmetric subspace
%  maxeig                - Maximum eigenvalue of a matrix    
%  mineig                - Minimum eigenvalue of a matrix 
%  trace2                - Trace-square of a matrix
%  trnorm                - Trace-norm
%  comm                  - Commutator
%  grstate               - Normalized ground state of a Hamiltonian
%  thstate               - Thermal ground state
%  thstate0              - Thermal ground state for T=0
%  addnoise              - Add white noise to a quantum state
%  binom                 - Binomial
%  qvec                  - Empty state vector for given number of qudits
%  qsize                 - Size of state vector or density matrix in qudits
%  qeye                  - Identity matrix for given number of qudits
%
%
% Commands with sparse matrices
%  spreordermat          - Sparse version of reordermat
%  spcoll                - Sparse version of coll
%  spinteract            - Sparse version of interact
%  spnnchain             - Sparse version of nnchain
%  spnnchainp            - Sparse version of nnchainp
%  spising               - Sparse version of ising
%  spisingp              - Sparse version of isningp
%  spquditop             - Sparse version of quditop
%  sptwoquditop          - Sparse version of twoquditop
%  splatticep            - Two-dimensional lattice Hamiltonian, periodic BC, sparse
%  splattice             - Two-dimensional lattice Hamiltonian, aperiodic BC, sparse
%  spising2Dp            - Two-dimensional Ising Hamiltonian, periodic BC, sparse
%
%
% Commands for symmetric states
%  sym2prodbasis         - Symmetric representation to product basis
%  sym2bipartite         - Symmetric representation to bipartite state
%  example_sym2bipartite - Example for using sym2bipartite
%  bipartite2prodbasis   - Bipartite representation to symmetric representation
%  ptsym                 - Partial transpose of a symmetric multiqubit state
%  example_ptsym         - Example for using ptsym
