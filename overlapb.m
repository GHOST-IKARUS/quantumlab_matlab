% overlapb   Maximal overlap with biseparable states
%            for a multi-qubit state. Use: overlapb(state)
%            where state is a state vector. The routine is
%            not based on numerical search thus gives
%            always a good result.

function maxoverlap=overlapb(state)

% Maximal overlap = the square of the maximum
% of the Schmidt coefficients corresponding to
% all bipartitioning

N=log2(length(state));

% Loop for all bipartitionings
maxoverlap=0;
for n=1:2^N-2
    list=[];
    b=dec2bin(n+2^N);
    for m=N:-1:1
       if b(m+1)=='1',
           list=[list,m];
       end %if
    end %for 
    if length(list)<=floor(N/2+0.51),
       maxoverlap=max(maxoverlap,max(eig(remove(ketbra(state),list))));
   end %if
end %for