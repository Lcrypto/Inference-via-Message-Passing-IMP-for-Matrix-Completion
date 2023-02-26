% A Code to test IMP algorithm on Netflix Data Matrix 2
clear all

% Select average number of observed ratings per user
n = 4; 

fprintf(1,'Generating Data Matrices...\n');
% Setup matrices and parameters
load data2

J = R>0; Je = J-Jt; 
[a b] = find(Je==1); 
tmp = sub2ind(size(Je),a,b); 
k = randperm(length(tmp)); 
n = n*size(R,1); indt = tmp(k(1:n));
Je = zeros(size(Je)); Je(indt) = 1; 

Re = Je.*R; Je = (Re>0); Rt = Jt.*R; 
gu = 16; gm = 16; 
clear R J

fprintf(1,'Learning with IMP...\n');
% Initialze via VDVQ clustering
[Pr] = imp1(Re,Je,gu,gm);   
% Do MP update
max_iter = 20;
[Pt Ps] = imp2(full(Re),Pr,gu,gm,max_iter);

fprintf(1,'Computing RMSE on Validataion Set...\n');
% Compute rmse on validation matrix
error = rmse(Rt,Jt,Pr,Ps,Pt,gu,gm);
fprintf(1,'RMSE: %f\n',error);