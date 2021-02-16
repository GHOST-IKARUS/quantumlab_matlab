% example_Fj_approx  Example to test Fj.m
%                    x axis: 0..1
%                    y axis: Fj(x) for j=14

dx=0.02;
xR=0.01:dx:0.99999;

j=28/2;

FjA=[];
Fj_approxA=[];
for x=xR
   FjA=[FjA,Fj(x,j)];
   Fj_approxA=[Fj_approxA,Fj_approx(x,j)];
end %if

plot(xR,FjA)
hold on
plot(xR,Fj_approxA,'--')
hold off
xlabel('x')
ylabel('F_j(x),F_{approx,j}(x)')
