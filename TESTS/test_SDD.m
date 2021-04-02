%-------------------------- test_SDD -------------------------------------%
%
% Script to test the "Shifted Diagonal Design (SDD)"
%
% In this example m=4, k=5
%
% The first 4 rows are expected to include
% (1,5,9,13);
% (2,6,10,14);
% (3,7,11,15);
% (4,8,12,16);
%-------------------------------------------------------------------------%
% 03/30/21, J.B., Initial version

clc;
clear;

% Adding paths to Algorithm and test matrix
addpath('../ALGS');

% Call SDD
pars.print=1;
m=4;
k=5;
[M,out]=SDD_PT(m,k,pars);

Mp1 = M(1:m,:);

% Test for m=3 (k=4)
m=3;
k=4;
[M3,out3]=SDD_PT(m,k,pars);