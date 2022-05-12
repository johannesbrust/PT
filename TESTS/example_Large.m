%-------------------------- example_Large --------------------------------%
%
% Script to test a "Finite Geometry Design (FGD)" and
% the COMP (Combinatorial orthogonal matching pursuit) algorithms
% to decode a set of pooled samples (e.g., corresponding to a COVID test)
%
% This problem is large with N = 100 000. A design based on
% pool sizes of 317 can be used to correctly identify upto 317 true
% positives.
%
% In this example 100 positives are identified
%
%-------------------------------------------------------------------------%
% 03/30/21, J.B., Initial version
% 05/12/22, J.B., Update to include the FGD algorithm

clc;
clear;

% Adding paths to Algorithm and test matrix
addpath('../ALGS');
addpath('../EXTERNAL');

% Call FGD (using an affine plane)
pars.print=1;
p = 317;
n = 1; % i.e., m=p^2=4
m = p^n;
k = 100;
[M,out]=FGD_PT(p,n,k,pars);

% Set defectives, 2,3,7,13
x = zeros(m*m,1);
x(1:1000:100000) = 1;

y = (sum(M(:,(x==1)),2)>0); %.*ones(t2,1);

% Call decoding algorithm
pars2.print = 1;
pars2.d = sum(x);

[x_c,out2] = COMP_PT((M>0),y,pars2);

err = norm(x-x_c);
