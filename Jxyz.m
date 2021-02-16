% Jxyz   Defines Jx, Jy, Jz collective angular momentum operators 
%        for N qubits. It can be used as [Jx,Jy,Jz]=Jxyz(N)

function [Jx,Jy,Jz]=Jxyz(N);

paulixyz;
Jx=coll(x,N)/2;
Jy=coll(y,N)/2;
Jz=coll(z,N)/2;