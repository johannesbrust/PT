%-------------------------- test_Pooling ---------------------------------%
%
% Script to test the "Shifted Diagonal Design (SDD)" and
% the COMP (Combinatorial orthogonal matching pursuit) algorithms
% to decode a set of pooled samples (e.g., corresponding to a COVID test)
%
% In this example m=4, k=3
%
%-------------------------------------------------------------------------%
% 03/30/21, J.B., Initial version

clc;
clear;

% Adding paths to Algorithm and test matrix
addpath('../ALGS');

% Call SDD
pars.print=1;
m=24;
k=7;
[M,out]=SDD_PT(m,k,pars);

% Set defectives, 2,3,13
x = zeros(m*m,1);
x(2)=1;
x(3)=1;
x(13)=1;
x(75)=1;
x(381)=1;

y = (sum(M(:,(x==1)),2)>0); %.*ones(t2,1);

% Call decoding algorithm
pars2.print = 1;
pars2.d = sum(x);

[x_c,out2] = COMP_PT((M>0),y,pars2);

err = norm(x-x_c);
