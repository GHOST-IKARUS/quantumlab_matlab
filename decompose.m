% decompose   Creates a string with the Pauli decomposition of a
%             Hermitian matrix. 
%    Form of use: decompose(rho) where rho is the density matrix.
%    One can print in LaTeX format with decompose(rho,1).
%    For the numbering of qubits we note that mkron(x,y,z) is
%    decomposed into "xyz" or "X^{(1)}Y^{(2)}Z^{(3)}".
%    For the LaTeX output the numbering is different from 
%    the one used in other commands (e.g., reorder) in which 
%    the numbering of qubits were 3 2 1 and not 1 2 3. This is 
%    to help to produce a LaTeX formula that can directly be
%    used in a paper. Form decompose(rho,mode,threshold) 
%    makes it possible to give a threshold. Below this 
%    value a correlation is considered zero. Default is 1e-14.

function pstring=decompose(rho,varargin)

% Correlations smaller than that are taken to be zero
mincorr=1e-14;

LaTeXON=0;
if ~isempty(varargin),
   if length(varargin)>=1,
      LaTeXON=varargin{1};
   end %if
   if length(varargin)==2,
      mincorr=varargin{2};
   end %if
   if length(varargin)>2,
      error('Wrong number of input arguments.')
   end %if
end %if
x=[0 1;1 0];
z=[1 0;0 -1];
y=i*x*z;
e=eye(2);

% Use sparse matrices
e=sparse(e);
x=sparse(x);
y=sparse(y);
z=sparse(z);

% rho=ketbra2(rho);
  
[sy,sx]=size(rho);

N=log2(sx);

pstring='';

OPstr='1xyz';
OPstr2='EXYZ';

opindex=zeros(N);

while 1
    switch opindex(1)
        case 0, op=e;
        case 1, op=x;
        case 2, op=y;
        case 3, op=z;
    end %switch 
    for m=2:N
        switch opindex(m)
            case 0, op=kron(op,e);
            case 1, op=kron(op,x);
            case 2, op=kron(op,y);
            case 3, op=kron(op,z);
        end %switch              
    end %for
    if trace(rho*op)~=0,  
        % Correlations: Use "real" to get rid off small imaginary parts
        corr=real(trace(rho*op))/2^N;
        termstr='';
        % Construct string for the term
        if LaTeXON,
            if opindex==zeros(N),
                termstr='E';
            else
                for m=1:N
                     if opindex (m) >0,                 
                        termstr=[termstr,OPstr2(opindex(m)+1) '^{(' num2str(m) ')}' ];  
                    end %if
                end %for
                
            end %if
        else
            for m=1:N
               termstr=[termstr,OPstr(opindex(m)+1)];
            end %for
        end %if
        if isempty(pstring),
            if abs(corr-1)<mincorr,  
                pstring=[termstr];
            elseif abs(corr+1)<mincorr,
                pstring=['-' termstr];
            elseif abs(corr)>mincorr,
                pstring=[num2str(corr) '*' termstr];
            end %if
        else
            if corr>mincorr,
                if abs(corr-1)<mincorr,
                    pstring=[pstring '+' termstr];
                else
                    pstring=[pstring '+' num2str(corr) '*' termstr];
                end %if
            elseif corr<-mincorr,
                if abs(corr+1)<mincorr,
                    pstring=[pstring '-' termstr];
                else
                    pstring=[pstring '-' num2str(abs(corr)) '*' termstr];
                end %if
            end %if             
        end %if
    end %if
    
    % Increase counter 
    k=N;
    opindex(k)=opindex(k)+1;
    while opindex(k)==4
        opindex(k)=0;
        k=k-1;
        if k==0,
            if isempty(pstring),
                pstring='0';
            end %if
            return
        end %if
        opindex(k)=opindex(k)+1;
    end %while
   
end %while

