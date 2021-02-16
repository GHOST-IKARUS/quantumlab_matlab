% Fj_approx   The function Fj(x) appearing in extreme spin squeezing 
%             of Sorensen and Molmer, using their approximate formula
%             See: http://arxiv.org/abs/quant-ph/0011035
%             Only for integer j and x<1.
%             Fj_approx(x,j) gives F_j(x). 
%             Fj_approx is a faster version of F_j, however, it gives
%             only a (very good) lower bound.
%           
%  See example_Fj_approx.m 

% Program by Geza Toth (C) 2016.

function Fjx=Fj_approx(x,j)

if x>1
    error('Implemented only for x<=1.');
end %if

% The expectation value of jx from the input
exJx_input=x*j;

% k-producabilty (k must be even!)
k=2*j;
%if ~even(k)
if mod(k,2)>0
    error('Implemented only for integer j.');
end %if

Fjx=((j+1)-j*x^2-sqrt(((j+1)-j*x^2)^2-x^2))/2;


