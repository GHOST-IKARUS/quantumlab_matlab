% Test new routines and compare them to old versions
% or test self-consistency

clear all
format compact

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% fisher.m
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Nit=100;
for n=1:Nit
    
    d=randi(10,1);
    r=rdmat(1,d);
    A=rdmat(1,d);
    B=rdmat(1,d);
    
    if norm(fisher(r,A)-fisher_V_5_5(r,A))>1e-14
        error('Error')
    end
    
    if norm(fisher(r,A,B)-fisher_V_5_5(r,A,B))>1e-14
        error('Error')
    end
    
end

% Not full rank
for n=1:Nit
    
    d=randi(10,1);
    d2=d-randi(d,1);
    r=blkdiag(rdmat(1,d2),zeros(d-d2,d-d2));
    A=rdmat(1,d);
    B=rdmat(1,d);
    
    if norm(fisher(r,A)-fisher_V_5_5(r,A))>1e-14
        error('Error')
    end
    
    if norm(fisher(r,A,B)-fisher_V_5_5(r,A,B))>1e-14
        error('Error')
    end
    
end

% Compare run times
d=4;
r=rdmat(1,d);
A=rdmat(1,d);
tic
for n=1:Nit   
    fisher(r,A);
end
toc_fisher=toc
tic
for n=1:Nit   
    fisher_V_5_5(r,A);
end
toc_fisher_V_5_5=toc
ratio_of_execution_times=toc_fisher_V_5_5/toc_fisher

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% pt.m
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
r=rdmat(4);

rpt1=pt(r,[4 3],2);
rpt2=pt(r,[2],4);
if norm(rpt1-rpt2)>1e-10
    error('Error') 
end

rpt1=pt_V_5_6(r,[4 3],2);
rpt2=pt(r,[3],[4 2 2]);
if norm(rpt1-rpt2)>1e-10
    error('Error') 
end

r=rdmat(8);
rpt1=pt(r,[8 7 6],2);
rpt2=pt(r,[2],[8 32]);
if norm(rpt1-rpt2)>1e-10
    error('Error') 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% reordermat
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

m1=reordermat_V_5_6([1 2],4);
m2=reordermat([2 1 3],[4 2 2]);
if norm(m1-m2)>1e-10
     error('Error') 
end

% DOES NOT WORK WELL
m1=reordermat_V_5_6([1 2 3],4);
m2=reordermat([2 1 4 3 5],[4 2 2 2 2]);
if norm(m1-m2)>1e-10
    error('Error') 
end

% I we do nothing then it is an idenity matrix
m1=reordermat([2 1],[2 4]);
m2=reordermat_V_5_6([3 2 1],2);
if norm(m1-m2)>1e-10
    error('Error') 
end

s1=rdmat(1,2);
s2=rdmat(1,3);
r1=mkron(s2,s1);
r2=reorder(r1,[2 1],[3 2]);
r2b=mkron(s2,s1);
if norm(r2-r2b)>1e-10
    error('Error') 
end

s1=rdmat(1,2);
s2=rdmat(1,3);
r1=mkron(s2,s1);
r2=reorder(r1,[1 2],[3 2]);
r2b=mkron(s1,s2);
if norm(r2-r2b)>1e-10
    error('Error') 
end

s1=rdmat(1,3);
s2=rdmat(1,2);
s3=rdmat(1,2);
r1=mkron(s3,s2,s1);
r2=reorder(r1,[2 3 1],[2 2 3]);
r2b=mkron(s2,s3,s1);
if norm(r2-r2b)>1e-10
    error('Error') 
end

s1=rdmat(1,5);
s2=rdmat(1,2);
s3=rdmat(1,3);
r1=mkron(s3,s2,s1);
r2=reorder(r1,[2 3 1],[3 2 5]); 
r2b=mkron(s2,s3,s1);
if norm(r2-r2b)>1e-10
    error('Error') 
end

s1=rdmat(1,2);
s2=rdmat(1,3);
s3=rdmat(1,2);
s4=rdmat(1,3);
r1=mkron(s4,s3,s2,s1);
r2=reorder(r1,[4 3 2 1],[2 2 3 3]); 
r2b=mkron(s4,s3,s2,s1);
if norm(r2-r2b)>1e-10
    error('Error') 
end

s1=rdmat(1,3);
s2=rdmat(1,2);
s3=rdmat(1,3);
s4=rdmat(1,2);
r1=mkron(s4,s3,s2,s1);
%r2=reorder(r1,[4 2 3 1],[2 2 3 3]); 
r2=reorder(r1,[4 2 3 1],[2 3 2 3]); 
r2b=mkron(s4,s2,s3,s1);
if norm(r2-r2b)>1e-10
    error('Error') 
end


