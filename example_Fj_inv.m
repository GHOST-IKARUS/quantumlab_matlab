% example_Fj:inv  Example to test Fj_inv.m
%                 x axis: 0..1
%                 y axis: Fj(x) for j=14

dx=0.02;
xR=0.01:dx:0.49;

j=28/2;

Fj_invA=[];
for x=xR
   Fj_invA=[Fj_invA,Fj_inv(x,j)];
end %if

plot(xR,Fj_invA)
xlabel('x')
ylabel('F_j^{-1}(x)')
