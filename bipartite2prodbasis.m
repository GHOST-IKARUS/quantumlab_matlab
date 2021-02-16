% bipartite2prodbasis    Converts a bipartite density matrix to
%                        a multiqubit state

function rho_prodbasis=sym2bipartite(rho_bipartite,M,N)

s1=M+1;
s2=N-M+1;

% Size of rho_bipartite
s=s1*s2;

rho_prodbasis=zeros(2^N,2^N);

for n=1:s
    n1=floor((n-1)/s2);
    n2=n-1-n1*s2;
    if n1>M || n2>N-M
        save
        error('Out of range.')
    end
    fn=kron(dstate(n1,M),dstate(n2,N-M));
    for m=1:s
        m1=floor((m-1)/s2);
        m2=m-1-m1*s2;
        if m1>M || m2>N-M
            save
            error('Out of range.')
        end
        fm=kron(dstate(m1,M),dstate(m2,N-M));
        rho_prodbasis=rho_prodbasis+rho_bipartite(n,m)*fn*fm.';
    end
end

