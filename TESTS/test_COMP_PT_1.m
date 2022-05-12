%-------------------- test_COMP_PT_1 -------------------------------------%
%
% Script to test Combinatorial Orthogonal Matching Pursuit (COMP) for
% Pooling Tests
%
% The test matrix is taken from Wikipedia
% https://en.wikipedia.org/wiki/Group_testing
% #Combinatorial_orthogonal_matching_pursuit_(COMP)
% 
% This is an example for which the algorithm identifies false positives
%-------------------------------------------------------------------------%
% 03/17/21, J.B., Initial version

clc;
clear;

% Adding paths to Algorithm and test matrix
addpath('../ALGS');
addpath('../AUXILIARY');

% Initialize RNG
rng(0);

% Generate matrix M
M = [0 1 1 0 0 0 0 0;
    1 0 0 0 0 0 0 0;
    0 0 0 1 0 1 0 1;
    0 0 1 0 0 0 0 0;
    0 0 0 1 0 0 0 0;
    0 1 0 0 0 0 1 1;
    1 0 0 1 0 0 0 0];

% Setup problem for M (and convert to binary)
M = (M>0);
[t,n] = size(M);

% Example x (with three positive cases)
x2 = zeros(n,1,'int8');
x2(2,1) = 1;
x2(6,1) = 1;
x2(7,1) = 1;

y2 = (sum(M(:,(x2==1)),2)>0);

% Call to algorithm
pars.print = 1;
pars.d = sum(x2);

[x2_c,out2] = COMP_PT(M,y2,pars);

% Call to DD algorithm
[x2_D,out2_D] = DD_PT(M,y2,pars);

% Call to SCOPM algorithm
[x2_S,out2_S] = SCOMP_PT(M,y2,pars);

