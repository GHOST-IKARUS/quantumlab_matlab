% example_Fj  Example to test Fj.m
%             x axis: 0..1
%             y axis: Fj(x) for j=14

dx=0.02;
xR=0.001:dx:0.999;

j=28/2;

FjA=[];
for x=xR
   FjA=[FjA,Fj(x,j)];
end %if

plot(xR,FjA)
xlabel('x')
ylabel('F_j(x)')
