% paulistr   Converts a Pauli string into an operator.
%   Form of use: paulistr(string). E.g., paulistr('5*xye+xyz') results
%   in 5*mkron(x,y,e)+mkron(x,y,z). Here x,y,z and e are defined by
%   the command paulixyz. The letter 'e' represents
%   the 2x2 identity matrix.

function op=paulistr(pstring)

x=[0 1;1 0];
z=[1 0;0 -1];
y=i*x*z;
e=eye(2);

index=1;

pstring2='';
while index<=length(pstring),
    s=pstring(index);
    if s=='e' || s=='x' || s=='y' || s=='z',
        pstring2=[pstring2,'mkron(' s];
        index=index+1;
        if index>length(pstring),
            pstring2=[pstring2 ')'];
	    op=eval(pstring2);
            return;
        end %if
        s=pstring(index);
        while s=='e' || s=='x' || s=='y' || s=='z',
           pstring2=[pstring2 ',' s];
           index=index+1;
           if index>length(pstring),
              pstring2=[pstring2 ')'];
	      op=eval(pstring2);
              return;
           end %if
           s=pstring(index);         
        end %while
        pstring2=[pstring2 ')'];
    else
        pstring2=[pstring2,s];
        index=index+1;
    end %if
end %while


